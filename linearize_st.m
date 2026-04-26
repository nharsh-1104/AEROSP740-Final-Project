function [A_d, B_d, c_d] = linearize_st(x_bar, u_bar, p, dt)
% LINEARIZE_ST  Numerical linearization of the RK4-discretized ST model.
%
%   [A_d, B_d, c_d] = linearize_st(x_bar, u_bar, p, dt)
%
% Computes the discrete affine model:
%   x_{k+1} ≈ A_d * x_k + B_d * u_k + c_d
%
% by numerically differentiating the RK4 integration step of the ST
% continuous dynamics.  This avoids error-prone hand-derivation of the
% 7-state Jacobian while retaining the accuracy of RK4.
%
% Inputs
%   x_bar  - Linearization state  [s_x; s_y; delta; v_x; psi; psi_dot; beta] (7x1)
%   u_bar  - Linearization input  [v_delta; a_x]                              (2x1)
%   p      - Vehicle parameter struct (parameters_vehicle2)
%   dt     - Sample time [s]
%
% Outputs
%   A_d   - Discrete state Jacobian   (7x7)
%   B_d   - Discrete input Jacobian   (7x2)
%   c_d   - Discrete affine offset    (7x1)

  nx = 7;
  nu = 2;

  % Nominal RK4 step
  x_next_nom = rk4_step_st(x_bar, u_bar, p, dt);

  % ── A_d: df_d / dx  via central finite differences ────────────────────
  %
  % Perturbation sizing is critical for the ST model.  States like psi_dot
  % and beta are often near zero, so a purely relative perturbation collapses
  % to machine-epsilon noise.  We use state-aware minimum perturbations:
  %
  %   s_x, s_y   : positions  — min 1e-4 m
  %   delta      : steering   — min 1e-4 rad  (~0.006 deg)
  %   v_x        : speed      — min 1e-3 m/s
  %   psi        : yaw        — min 1e-4 rad
  %   psi_dot    : yaw rate   — min 1e-3 rad/s
  %   beta       : side-slip  — min 1e-3 rad
  %
  eps_floor = [1e-4; 1e-4; 1e-4; 1e-3; 1e-4; 1e-3; 1e-3];
  eps_x = max(eps_floor, 1e-4 * abs(x_bar));

  A_d = zeros(nx, nx);
  for j = 1:nx
    x_plus  = x_bar;  x_plus(j)  = x_plus(j)  + eps_x(j);
    x_minus = x_bar;  x_minus(j) = x_minus(j) - eps_x(j);

    f_plus  = rk4_step_st(x_plus,  u_bar, p, dt);
    f_minus = rk4_step_st(x_minus, u_bar, p, dt);

    A_d(:, j) = (f_plus - f_minus) / (2 * eps_x(j));
  end

  % ── B_d: df_d / du  via central finite differences ────────────────────
  %   v_delta : min 1e-3 rad/s
  %   a_x     : min 1e-3 m/s^2
  B_d = zeros(nx, nu);
  eps_u = max(1e-3 * ones(nu,1), 1e-4 * abs(u_bar));

  for j = 1:nu
    u_plus  = u_bar;  u_plus(j)  = u_plus(j)  + eps_u(j);
    u_minus = u_bar;  u_minus(j) = u_minus(j) - eps_u(j);

    f_plus  = rk4_step_st(x_bar, u_plus,  p, dt);
    f_minus = rk4_step_st(x_bar, u_minus, p, dt);

    B_d(:, j) = (f_plus - f_minus) / (2 * eps_u(j));
  end

  % ── Affine offset ─────────────────────────────────────────────────────
  % Ensures A_d*x_bar + B_d*u_bar + c_d = f_d(x_bar, u_bar) exactly.
  c_d = x_next_nom - A_d * x_bar - B_d * u_bar;
end


% ═══════════════════════════════════════════════════════════════════════════
function x_next = rk4_step_st(x, u, p, dt)
% RK4_STEP_ST  Single RK4 integration step of ST continuous dynamics.
  k1 = st_dynamics_continuous(x,              u, p);
  k2 = st_dynamics_continuous(x + dt/2 * k1,  u, p);
  k3 = st_dynamics_continuous(x + dt/2 * k2,  u, p);
  k4 = st_dynamics_continuous(x + dt   * k3,  u, p);
  x_next = x + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
end
