function res = run_nonlinear_mpc(x_ks0, x_mb0, N_sim, hp, hu, mpc, dlc, p, v0)
% RUN_NONLINEAR_MPC  Nonlinear MPC (single-shooting) with KS predictor + MB plant.
%
%   res = run_nonlinear_mpc(x_ks0, x_mb0, N_sim, hp, hu, mpc, dlc, p, v0)
%
% Architecture
%   Predictor:  CommonRoad Kinematic Single-Track (KS) model, Euler-discretised
%               x = [s_x, s_y, delta, v, psi];   u = [v_delta, a_x]
%   Plant:      CommonRoad Multi-Body (MB) model — vehicleDynamics_MB
%               29 states, Pacejka Magic Formula combined-slip tires
%   Formulation: Single-shooting NLP (ONLY control sequence U is optimised).
%               Initial state enters as a parameter, states are propagated
%               symbolically inside CasADi — no equality constraints needed.
%   Solver:     IPOPT via CasADi
%
% The NLP is assembled ONCE; at each of the N_sim simulation steps only
% the parameter vector (current state + reference) is updated.
%
% Inputs
%   x_ks0  - Initial KS state  [s_x; s_y; delta; v; psi]  (5×1)
%   x_mb0  - Initial MB state  (29×1, from init_MB)
%   N_sim  - Number of simulation steps
%   hp     - Prediction horizon
%   hu     - Control horizon  (hu ≤ hp)
%   mpc    - Struct: .dt, .nx, .nu, .u_min, .u_max, .d_min, .d_max,
%                    .Q, .R, .QN
%   dlc    - DLC path parameter struct
%   p      - Vehicle2 parameter struct (from parameters_vehicle2)
%   v0     - Reference longitudinal speed [m/s]
%
% Output  res
%   .sx, .sy, .psi, .v  - State trajectories  (1 × N_sim+1)
%   .u                  - Applied control inputs  (2 × N_sim)
%   .tc                 - Per-step solve times  (1 × N_sim)
%   .t_avg, .t_max      - Timing statistics [s]
%   .sy_ref             - Reference sy at each recorded sx  (1 × N_sim+1)
%   .eps                - Lateral RMSE [m]

  import casadi.*;

  dt = mpc.dt;
  nx = mpc.nx;   % 5
  nu = mpc.nu;   % 2
  l  = mpc.l;    % wheelbase from paper (2.5 m)

  % ──────────────────────────────────────────────────────────────────────
  % Build CasADi NLP (single shooting, done once)
  %
  % Decision variables: U = [u_0; u_1; ...; u_{hu-1}]   (nu*hu × 1)
  %
  % Parameters:   P = [x0  (nx);   Xref  (nx × (hp+1) stacked)]
  %   Xref(:,1)  = reference at step 0 (current position, not used in cost)
  %   Xref(:,k+1) = reference at step k  for k = 1..hp
  %
  % Steering-angle state constraints: delta_min ≤ delta_k ≤ delta_max
  %   collected as g_delta  (hp × 1)
  % ──────────────────────────────────────────────────────────────────────

  n_Uvars = nu * hu;
  n_Pvars = nx * (hp + 2);   % x0 (nx) + Xref (nx*(hp+1))

  U_sym = MX.sym('U', n_Uvars);
  P_sym = MX.sym('P', n_Pvars);

  x0_p   = P_sym(1:nx);
  Xref_p = reshape(P_sym(nx + 1 : end), nx, hp + 1);  % (5 × hp+1)

  cost       = MX(0);
  g_delta    = MX(hp, 1);    % predicted steering angle at each step

  x_k = x0_p;    % symbolic propagation starts from initial state

  for k = 1:hp
    % Select the appropriate control (hold last after control horizon)
    u_start = (min(k, hu) - 1) * nu + 1;
    u_k     = U_sym(u_start : u_start + nu - 1);

    % Euler-discretised KS dynamics
    x_k = x_k + dt * [ x_k(4)*cos(x_k(5));
                        x_k(4)*sin(x_k(5));
                        u_k(1);
                        u_k(2);
                        x_k(4)/l * tan(x_k(3)) ];

    % Collect predicted steering angle for constraint
    g_delta(k) = x_k(3);

    % Cost (stage or terminal)
    e = x_k - Xref_p(:, k + 1);
    if k < hp
      cost = cost + e' * mpc.Q * e;
    else
      cost = cost + e' * mpc.QN * e;
    end
    cost = cost + u_k' * mpc.R * u_k;
  end

  % Solver bounds
  lbx = repmat(mpc.u_min, hu, 1);
  ubx = repmat(mpc.u_max, hu, 1);
  lbg = mpc.d_min * ones(hp, 1);   % steering angle lower bound at each step
  ubg = mpc.d_max * ones(hp, 1);   % steering angle upper bound

  ipopt_opts                    = struct();
  ipopt_opts.ipopt.print_level  = 0;
  ipopt_opts.ipopt.max_iter     = 500;
  ipopt_opts.ipopt.tol          = 1e-6;
  ipopt_opts.print_time         = 0;

  solver = nlpsol('nmpc_solver', 'ipopt', ...
                  struct('x', U_sym, 'f', cost, 'g', g_delta, 'p', P_sym), ...
                  ipopt_opts);

  % ── Simulation loop ───────────────────────────────────────────────────
  sx_h  = zeros(1, N_sim + 1);
  sy_h  = zeros(1, N_sim + 1);
  psi_h = zeros(1, N_sim + 1);
  v_h   = zeros(1, N_sim + 1);
  u_h   = zeros(nu, N_sim);
  tc_h  = zeros(1, N_sim);

  x_ks = x_ks0(:);
  x_mb = x_mb0(:);

  sx_h(1)  = x_ks(1);   sy_h(1)  = x_ks(2);
  psi_h(1) = x_ks(5);   v_h(1)   = x_ks(4);

  U0 = zeros(n_Uvars, 1);   % warm-start (updated each step)

  for k = 1:N_sim
    % Predict future sx (constant speed assumption)
    sx_fut = x_ks(1) + v0 * (0:hp) * dt;

    % Reference over horizon
    [sy_ref_v, psi_ref_v] = get_dlc_reference(sx_fut, dlc);

    Xref_num = zeros(nx, hp + 1);
    for j = 1:hp + 1
      Xref_num(:, j) = [sx_fut(j); sy_ref_v(j); 0; v0; psi_ref_v(j)];
    end

    P_num = [x_ks; Xref_num(:)];

    % Solve NLP
    t_s = tic;
    sol = solver('x0',  U0,    'p',   P_num, ...
                 'lbx', lbx,   'ubx', ubx, ...
                 'lbg', lbg,   'ubg', ubg);
    tc_h(k) = toc(t_s);

    U_opt = full(sol.x);
    u_k   = U_opt(1:nu);

    % Shift warm-start: drop first input, append zero
    U0 = [U_opt(nu + 1 : end); zeros(nu, 1)];

    % Apply to MB plant
    x_mb = simulate_plant_mb(x_mb, u_k, p, dt);
    x_ks = mb_to_ks(x_mb);

    % Store
    sx_h(k+1)  = x_ks(1);   sy_h(k+1)  = x_ks(2);
    psi_h(k+1) = x_ks(5);   v_h(k+1)   = x_ks(4);
    u_h(:, k)  = u_k;
  end

  % ── Performance metrics ───────────────────────────────────────────────
  [sy_ref_all, ~] = get_dlc_reference(sx_h, dlc);

  % Reference sx trajectory (constant speed assumption for reference)
  sx_ref_all = sx_h(1) + v0 * (0:N_sim) * dt;

  res.sx     = sx_h;
  res.sy     = sy_h;
  res.sx_ref = sx_ref_all;
  res.sy_ref = sy_ref_all;
  res.psi    = psi_h;
  res.v      = v_h;
  res.u      = u_h;
  res.tc     = tc_h;
  res.t_avg  = mean(tc_h);
  res.t_max  = max(tc_h);
  % Paper eq. (24): eps = (1/2n) * sum( (px-px_ref)^2 + (py-py_ref)^2 )
  n_pts   = numel(sx_h);
  res.eps = (1/(2*n_pts)) * sum((sx_h - sx_ref_all).^2 + (sy_h - sy_ref_all).^2);
end
