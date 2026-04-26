function res = run_nonlinear_mpc_st(x_st0, x_mb0, N_sim, hp, hu, mpc, dlc, p, v0)
% RUN_NONLINEAR_MPC_ST  Nonlinear MPC with ST predictor (RK4) + MB plant.
%
%   res = run_nonlinear_mpc_st(x_st0, x_mb0, N_sim, hp, hu, mpc, dlc, p, v0)
%
% Architecture
%   Predictor:  CommonRoad Single-Track (ST) dynamic bicycle model (7 states)
%               x = [s_x, s_y, delta, v_x, psi, psi_dot, beta]
%               u = [v_delta, a_x]
%               Includes lateral tire forces via linearized Pacejka cornering
%               stiffness, load transfer, yaw & side-slip dynamics.
%   Integration: 4th-order Runge-Kutta (RK4) — symbolic, inside CasADi
%   Plant:      CommonRoad Multi-Body (MB) model (29 states)
%   Formulation: Single-shooting NLP, solved by IPOPT
%
% The NLP is assembled ONCE; at each simulation step only the parameter
% vector (current state + reference) is updated.
%
% Inputs
%   x_st0  - Initial ST state  (7x1)
%   x_mb0  - Initial MB state  (29x1)
%   N_sim  - Number of simulation steps
%   hp     - Prediction horizon
%   hu     - Control horizon (hu <= hp)
%   mpc    - Struct with .dt, .nx_st, .nu, .u_min, .u_max, .d_min, .d_max,
%            .Q_st, .R, .QN_st
%   dlc    - DLC path parameter struct
%   p      - Vehicle parameter struct
%   v0     - Reference longitudinal speed [m/s]
%
% Output  res  (same fields as KS-based version for easy comparison)

  import casadi.*;

  dt  = mpc.dt;
  nx  = mpc.nx_st;   % 7
  nu  = mpc.nu;      % 2

  % ── Vehicle parameters needed for ST dynamics ─────────────────────────
  g   = 9.81;
  mu  = p.tire.p_dy1;
  C_Sf = -p.tire.p_ky1 / p.tire.p_dy1;
  C_Sr = -p.tire.p_ky1 / p.tire.p_dy1;
  lf  = p.a;
  lr  = p.b;
  lwb = lf + lr;
  h_s = p.h_s;
  m   = p.m;
  I_z = p.I_z;

  % ──────────────────────────────────────────────────────────────────────
  %  CasADi symbolic ST dynamics (for v_x > 0.5)
  % ──────────────────────────────────────────────────────────────────────
  x_sym = MX.sym('x', nx);
  u_sym = MX.sym('u', nu);

  % Unpack
  s_x_    = x_sym(1);  s_y_    = x_sym(2);
  delta_  = x_sym(3);  v_x_    = x_sym(4);
  psi_    = x_sym(5);  psi_d_  = x_sym(6);
  beta_   = x_sym(7);
  v_delta_ = u_sym(1); a_x_    = u_sym(2);

  % Gravity + load-transfer terms (NOT actual forces — matches CommonRoad
  % vehicleDynamics_ST which uses these raw terms directly).
  glr_ = g * lr - a_x_ * h_s;    % front-axle related  [m^2/s^2]
  glf_ = g * lf + a_x_ * h_s;    % rear-axle related   [m^2/s^2]

  % Continuous dynamics  (safe-guarded v_x for CasADi)
  v_x_safe = if_else(v_x_ > 0.5, v_x_, 0.5);

  f_sym = vertcat( ...
    v_x_ * cos(beta_ + psi_), ...                                    % dot s_x
    v_x_ * sin(beta_ + psi_), ...                                    % dot s_y
    v_delta_, ...                                                      % dot delta
    a_x_, ...                                                          % dot v_x
    psi_d_, ...                                                        % dot psi
    -mu*m/(v_x_safe*I_z*lwb)*(lf^2*C_Sf*glr_ + lr^2*C_Sr*glf_)*psi_d_ ...
      + mu*m/(I_z*lwb)*(lr*C_Sr*glf_ - lf*C_Sf*glr_)*beta_ ...
      + mu*m/(I_z*lwb)*lf*C_Sf*glr_*delta_, ...                        % dot psi_dot
    (mu/(v_x_safe^2*lwb)*(C_Sr*glf_*lr - C_Sf*glr_*lf) - 1)*psi_d_ ...
      - mu/(v_x_safe*lwb)*(C_Sr*glf_ + C_Sf*glr_)*beta_ ...
      + mu/(v_x_safe*lwb)*C_Sf*glr_*delta_ );                          % dot beta

  f_fun = Function('f_st', {x_sym, u_sym}, {f_sym});

  % ── RK4 integrator step (symbolic) ────────────────────────────────────
  k1_ = f_fun(x_sym, u_sym);
  k2_ = f_fun(x_sym + dt/2 * k1_, u_sym);
  k3_ = f_fun(x_sym + dt/2 * k2_, u_sym);
  k4_ = f_fun(x_sym + dt   * k3_, u_sym);
  x_next_sym = x_sym + (dt/6) * (k1_ + 2*k2_ + 2*k3_ + k4_);

  rk4_step = Function('rk4_st', {x_sym, u_sym}, {x_next_sym});

  % ──────────────────────────────────────────────────────────────────────
  %  Build NLP (single shooting)
  % ──────────────────────────────────────────────────────────────────────
  n_Uvars = nu * hu;
  n_Pvars = nx * (hp + 2);   % x0 (nx) + Xref (nx*(hp+1))

  U_opt_sym = MX.sym('U', n_Uvars);
  P_sym     = MX.sym('P', n_Pvars);

  x0_p   = P_sym(1:nx);
  Xref_p = reshape(P_sym(nx+1:end), nx, hp+1);   % (7 x hp+1)

  cost    = MX(0);
  g_delta = MX(hp, 1);
  x_k     = x0_p;

  for k = 1:hp
    % Select control (hold last after hu)
    u_idx = (min(k, hu) - 1) * nu + 1;
    u_k   = U_opt_sym(u_idx : u_idx + nu - 1);

    % RK4 propagation
    x_k = rk4_step(x_k, u_k);

    % Steering angle constraint
    g_delta(k) = x_k(3);

    % Tracking error cost
    e = x_k - Xref_p(:, k+1);
    if k < hp
      cost = cost + e' * mpc.Q_st * e;
    else
      cost = cost + e' * mpc.QN_st * e;
    end
    cost = cost + u_k' * mpc.R * u_k;
  end

  % Solver
  lbx = repmat(mpc.u_min, hu, 1);
  ubx = repmat(mpc.u_max, hu, 1);
  lbg = mpc.d_min * ones(hp, 1);
  ubg = mpc.d_max * ones(hp, 1);

  opts                    = struct();
  opts.ipopt.print_level  = 0;
  opts.ipopt.max_iter     = 500;
  opts.ipopt.tol          = 1e-6;
  opts.print_time         = 0;

  solver = nlpsol('nmpc_st', 'ipopt', ...
                  struct('x', U_opt_sym, 'f', cost, 'g', g_delta, 'p', P_sym), ...
                  opts);

  % ── Simulation loop ───────────────────────────────────────────────────
  sx_h  = zeros(1, N_sim+1);
  sy_h  = zeros(1, N_sim+1);
  psi_h = zeros(1, N_sim+1);
  v_h   = zeros(1, N_sim+1);
  u_h   = zeros(nu, N_sim);
  tc_h  = zeros(1, N_sim);

  x_st = x_st0(:);
  x_mb = x_mb0(:);

  sx_h(1) = x_st(1);  sy_h(1)  = x_st(2);
  psi_h(1)= x_st(5);  v_h(1)   = x_st(4);

  U0 = zeros(n_Uvars, 1);

  for k = 1:N_sim
    % Future sx (constant speed assumption for reference lookup)
    sx_fut = x_st(1) + v0 * (0:hp) * dt;

    % DLC reference
    [sy_ref_v, psi_ref_v] = get_dlc_reference(sx_fut, dlc);

    % Curvature-based steady-state references for psi_dot and beta
    [psid_ref_v, beta_ref_v] = get_curvature_references(sx_fut, dlc, v0, p.b);

    % Build reference for 7-state model:
    %   [sx_ref, sy_ref, 0(delta), v0, psi_ref, psi_dot_ref, beta_ref]
    Xref_num = zeros(nx, hp+1);
    for j = 1:hp+1
      Xref_num(:,j) = [sx_fut(j); sy_ref_v(j); 0; v0; psi_ref_v(j); psid_ref_v(j); beta_ref_v(j)];
    end

    P_num = [x_st; Xref_num(:)];

    % Solve NLP
    t_s = tic;
    sol = solver('x0', U0, 'p', P_num, ...
                 'lbx', lbx, 'ubx', ubx, ...
                 'lbg', lbg, 'ubg', ubg);
    tc_h(k) = toc(t_s);

    U_star = full(sol.x);
    u_k    = U_star(1:nu);

    % Warm-start shift
    U0 = [U_star(nu+1:end); zeros(nu,1)];

    % Apply to MB plant
    x_mb = simulate_plant_mb(x_mb, u_k, p, dt);
    x_st = mb_to_st(x_mb);

    % Store
    sx_h(k+1) = x_st(1);   sy_h(k+1) = x_st(2);
    psi_h(k+1)= x_st(5);   v_h(k+1)  = x_st(4);
    u_h(:,k)  = u_k;
  end

  % ── Performance metrics ───────────────────────────────────────────────
  [sy_ref_all, ~] = get_dlc_reference(sx_h, dlc);
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
  n_pts   = numel(sx_h);
  res.eps = (1/(2*n_pts)) * sum((sx_h - sx_ref_all).^2 + (sy_h - sy_ref_all).^2);
end
