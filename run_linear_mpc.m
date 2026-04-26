function res = run_linear_mpc(x_ks0, x_mb0, N_sim, hp, hu, mpc, dlc, p, v0)
% RUN_LINEAR_MPC  Linear time-varying MPC using KS predictor and MB plant.
%
%   res = run_linear_mpc(x_ks0, x_mb0, N_sim, hp, hu, mpc, dlc, p, v0)
%
% Architecture
%   Predictor:  KS model, first-order Taylor linearized at current state
%               then Euler-discretized → affine model per step
%   Plant:      CommonRoad Multi-Body (MB) model (vehicleDynamics_MB)
%   Solver:     Condensed QP solved with MATLAB's quadprog
%               (Optimization Toolbox required)
%
% At each time step:
%   1. Linearize KS model at (x_ks, u_bar=0) → A_d, B_d, c_d
%   2. Build condensed prediction matrices Phi, Gamma, D
%   3. Formulate and solve QP for optimal control sequence U*
%   4. Apply first control u*_0 to MB plant
%   5. Observe [sx, sy, delta, vx, psi] from MB state
%
% Cost functional (condensed form):
%   J = (X_pred − X_ref)' Q_bar (X_pred − X_ref)  +  U' R_bar U
%
% Inputs / Outputs  — identical signature to run_nonlinear_mpc.m

  dt = mpc.dt;
  nx = mpc.nx;   % 5
  nu = mpc.nu;   % 2
  l  = mpc.l;    % wheelbase from paper (2.5 m)

  % ── Block-diagonal cost matrices ─────────────────────────────────────
  % Stage weight on steps 1..(hp-1), terminal weight on step hp
  if hp == 1
    Q_bar = mpc.QN;
  else
    Q_bar = blkdiag(kron(eye(hp - 1), mpc.Q), mpc.QN);  % (hp*nx × hp*nx)
  end

  % Paper eq. (5): J = sum_{i=0}^{N-1} (x_i'Qx_i + u_i'Ru_i) + x_N'S_Nx_N
  % When hu < hp, the held last input u_{hu-1} is applied at steps hu..hp,
  % so it must be penalised (hp - hu + 1) times, not just once.
  % Build R_bar via the reduction matrix S that maps hu decision vars to hp inputs.
  S_red = zeros(hp * nu, hu * nu);
  for j = 1:hp
    j_hu = min(j, hu);
    S_red((j-1)*nu+(1:nu), (j_hu-1)*nu+(1:nu)) = eye(nu);
  end
  R_bar = S_red' * kron(eye(hp), mpc.R) * S_red;         % (hu*nu × hu*nu)

  % Symmetrize (numerical safety)
  Q_bar = (Q_bar + Q_bar') / 2;
  R_bar = (R_bar + R_bar') / 2;

  % quadprog options
  qp_opts = optimoptions('quadprog', ...
                         'Display',            'off', ...
                         'MaxIterations',       500,  ...
                         'OptimalityTolerance', 1e-8);

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

  lb = repmat(mpc.u_min, hu, 1);   % input lower bounds (hu*nu × 1)
  ub = repmat(mpc.u_max, hu, 1);   % input upper bounds

  u_prev = zeros(nu, 1);           % warm-start for quadprog

  for k = 1:N_sim
    % ── Reference trajectory ─────────────────────────────────────────
    sx_fut = x_ks(1) + v0 * (1:hp) * dt;   % predict ahead (steps 1..hp)
    [sy_ref_v, psi_ref_v] = get_dlc_reference(sx_fut, dlc);

    % Stack as X_ref (hp*nx × 1), reference for predicted steps 1..hp
    X_ref = zeros(hp * nx, 1);
    for j = 1:hp
      X_ref((j-1)*nx + (1:nx)) = ...
          [sx_fut(j); sy_ref_v(j); 0; v0; psi_ref_v(j)];
    end

    % ── Linearize KS model at current state ──────────────────────────
    u_bar = zeros(nu, 1);    % linearize at rest input
    [A_d, B_d, c_d] = linearize_ks(x_ks, u_bar, l, dt);

    % ── Build condensed prediction matrices ──────────────────────────
    % X_pred = Phi*x_ks + Gamma*U + D
    [Phi, Gamma, D] = build_prediction_matrices(A_d, B_d, c_d, nx, nu, hp, hu);

    % ── QP matrices ──────────────────────────────────────────────────
    %   min  0.5*U'*H_qp*U + f_qp'*U
    bias  = Phi * x_ks + D - X_ref;   % (hp*nx × 1)
    H_qp  = 2 * (Gamma' * Q_bar * Gamma + R_bar);
    f_qp  = 2 * Gamma' * Q_bar * bias;

    % Symmetrize H for numerical robustness
    H_qp = (H_qp + H_qp') / 2;

    % ── Steering angle inequality constraints ────────────────────────
    % Predicted delta at step j:  Phi(row_j,:)*x0 + Gamma(row_j,:)*U + D(row_j)
    % Extract rows corresponding to delta (state index 3, 1-based)
    delta_rows = 3 : nx : hp * nx;          % indices of delta in X_pred

    C_delta    = Gamma(delta_rows, :);                    % (hp × hu*nu)
    delta_bias = Phi(delta_rows, :) * x_ks + D(delta_rows);  % (hp × 1)

    % d_min ≤ C_delta*U + delta_bias ≤ d_max  →  2-sided inequality
    A_ineq = [ C_delta; -C_delta];
    b_ineq = [ mpc.d_max * ones(hp, 1) - delta_bias; ...
               delta_bias - mpc.d_min * ones(hp, 1) ];

    % ── Solve QP ─────────────────────────────────────────────────────
    U0_warm = repmat(u_prev, hu, 1);   % warm-start

    t_s = tic;
    [U_opt, ~, flag] = quadprog(H_qp, f_qp, A_ineq, b_ineq, ...
                                [], [], lb, ub, U0_warm, qp_opts);
    tc_h(k) = toc(t_s);

    if isempty(U_opt) || flag < 0
      warning('LMPC step %d: quadprog infeasible (flag=%d). Using previous input.', ...
              k, flag);
      U_opt = U0_warm;
    end

    u_k    = U_opt(1:nu);
    u_prev = u_k;

    % ── Apply to MB plant ─────────────────────────────────────────────
    x_mb = simulate_plant_mb(x_mb, u_k, p, dt);
    x_ks = mb_to_ks(x_mb);

    % ── Store ─────────────────────────────────────────────────────────
    sx_h(k+1)  = x_ks(1);   sy_h(k+1)  = x_ks(2);
    psi_h(k+1) = x_ks(5);   v_h(k+1)   = x_ks(4);
    u_h(:, k)  = u_k;
  end

  % ── Metrics ──────────────────────────────────────────────────────────
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
