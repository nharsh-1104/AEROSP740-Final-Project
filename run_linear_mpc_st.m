function res = run_linear_mpc_st(x_st0, x_mb0, N_sim, hp, hu, mpc, dlc, p, v0)
% RUN_LINEAR_MPC_ST  Linear time-varying MPC with ST predictor + MB plant.
%
%   res = run_linear_mpc_st(x_st0, x_mb0, N_sim, hp, hu, mpc, dlc, p, v0)
%
% Architecture
%   Predictor:  ST model, numerically linearized via RK4-based Jacobian
%               at current state each step -> A_d, B_d, c_d
%   Plant:      CommonRoad Multi-Body (MB) model (29 states)
%   Solver:     Condensed QP solved with MATLAB's quadprog
%
% State:  x = [s_x, s_y, delta, v_x, psi, psi_dot, beta]   (7 x 1)
% Input:  u = [v_delta, a_x]                                 (2 x 1)
%
% Inputs / Outputs — same structure as run_nonlinear_mpc_st.m

  dt = mpc.dt;
  nx = mpc.nx_st;   % 7
  nu = mpc.nu;      % 2

  % ── Block-diagonal cost matrices ─────────────────────────────────────
  if hp == 1
    Q_bar = mpc.QN_st;
  else
    Q_bar = blkdiag(kron(eye(hp-1), mpc.Q_st), mpc.QN_st);
  end

  % Input cost: penalise held inputs correctly (same as LMPC KS fix)
  S_red = zeros(hp * nu, hu * nu);
  for j = 1:hp
    j_hu = min(j, hu);
    S_red((j-1)*nu+(1:nu), (j_hu-1)*nu+(1:nu)) = eye(nu);
  end
  R_bar = S_red' * kron(eye(hp), mpc.R) * S_red;

  Q_bar = (Q_bar + Q_bar') / 2;
  R_bar = (R_bar + R_bar') / 2;

  % quadprog options
  qp_opts = optimoptions('quadprog', ...
                         'Display',            'off', ...
                         'MaxIterations',       500,  ...
                         'OptimalityTolerance', 1e-8);

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

  lb = repmat(mpc.u_min, hu, 1);
  ub = repmat(mpc.u_max, hu, 1);
  u_prev = zeros(nu, 1);

  for k = 1:N_sim
    % ── Reference trajectory ─────────────────────────────────────────
    sx_fut = x_st(1) + v0 * (1:hp) * dt;
    [sy_ref_v, psi_ref_v] = get_dlc_reference(sx_fut, dlc);

    % Curvature-based steady-state references for psi_dot and beta
    [psid_ref_v, beta_ref_v] = get_curvature_references(sx_fut, dlc, v0, p.b);

    X_ref = zeros(hp * nx, 1);
    for j = 1:hp
      X_ref((j-1)*nx + (1:nx)) = ...
          [sx_fut(j); sy_ref_v(j); 0; v0; psi_ref_v(j); psid_ref_v(j); beta_ref_v(j)];
    end

    % ── Linearize ST model at current state ──────────────────────────
    u_bar = zeros(nu, 1);
    [A_d, B_d, c_d] = linearize_st(x_st, u_bar, p, dt);

    % ── Build condensed prediction matrices ──────────────────────────
    [Phi, Gamma, D] = build_prediction_matrices(A_d, B_d, c_d, nx, nu, hp, hu);

    % ── QP matrices ──────────────────────────────────────────────────
    bias  = Phi * x_st + D - X_ref;
    H_qp  = 2 * (Gamma' * Q_bar * Gamma + R_bar);
    f_qp  = 2 * Gamma' * Q_bar * bias;
    H_qp  = (H_qp + H_qp') / 2;

    % Small Tikhonov regularization for numerical robustness
    H_qp  = H_qp + 1e-6 * eye(size(H_qp));

    % ── Steering angle inequality constraints ────────────────────────
    delta_rows = 3 : nx : hp * nx;
    C_delta    = Gamma(delta_rows, :);
    delta_bias = Phi(delta_rows, :) * x_st + D(delta_rows);

    A_ineq = [ C_delta; -C_delta];
    b_ineq = [ mpc.d_max * ones(hp,1) - delta_bias; ...
               delta_bias - mpc.d_min * ones(hp,1) ];

    % ── Solve QP ─────────────────────────────────────────────────────
    U0_warm = repmat(u_prev, hu, 1);

    t_s = tic;
    [U_opt, ~, flag] = quadprog(H_qp, f_qp, A_ineq, b_ineq, ...
                                [], [], lb, ub, U0_warm, qp_opts);
    tc_h(k) = toc(t_s);

    if isempty(U_opt) || flag < 0
      warning('LMPC_ST step %d: quadprog flag=%d. Using previous input.', k, flag);
      U_opt = U0_warm;
    end

    u_k    = U_opt(1:nu);
    u_prev = u_k;

    % ── Apply to MB plant ─────────────────────────────────────────────
    x_mb = simulate_plant_mb(x_mb, u_k, p, dt);
    x_st = mb_to_st(x_mb);

    % ── Store ─────────────────────────────────────────────────────────
    sx_h(k+1) = x_st(1);   sy_h(k+1) = x_st(2);
    psi_h(k+1)= x_st(5);   v_h(k+1)  = x_st(4);
    u_h(:,k)  = u_k;
  end

  % ── Metrics ──────────────────────────────────────────────────────────
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
