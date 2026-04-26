function [A_d, B_d, c_d] = linearize_ks(x_bar, u_bar, l, dt)
% LINEARIZE_KS  First-order Taylor linearization of the KS bicycle model,
%               discretized via the Euler method.
%
%   [A_d, B_d, c_d] = linearize_ks(x_bar, u_bar, l, dt)
%
% Continuous KS dynamics (CommonRoad vehicleDynamics_KS, no constraint
% clamping — bounds are handled by the MPC optimizer):
%
%   dot_s_x   = v * cos(psi)
%   dot_s_y   = v * sin(psi)
%   dot_delta = v_delta              (u_1)
%   dot_v     = a_x                  (u_2)
%   dot_psi   = v / l * tan(delta)
%
% Linearized around (x_bar, u_bar) and Euler-discretized:
%
%   x_{k+1} ≈  A_d * x_k  +  B_d * u_k  +  c_d
%
% where the affine offset c_d ensures consistency at the linearization
% point (i.e. the Euler step of the nonlinear model is reproduced exactly
% when x = x_bar, u = u_bar).
%
% Inputs
%   x_bar  - Linearization state  [s_x; s_y; delta; v; psi]  (5x1)
%   u_bar  - Linearization input  [v_delta; a_x]              (2x1)
%   l      - Wheelbase [m]
%   dt     - Sample time [s]
%
% Outputs
%   A_d   - Discrete state matrix   (5x5)
%   B_d   - Discrete input matrix   (5x2)
%   c_d   - Discrete affine offset  (5x1)

  delta_b = x_bar(3);
  v_b     = x_bar(4);
  psi_b   = x_bar(5);

  % Guard against near-singular steering angle (|delta| ≈ pi/2 is unphysical
  % for a car but add a tiny regulariser for robustness).
  cos_d2 = max(cos(delta_b)^2, 1e-6);

  % ── Continuous Jacobians ─────────────────────────────────────────────
  % A_c = df/dx  evaluated at (x_bar, u_bar)
  %         s_x   s_y   delta              v                psi
  Ac = [ 0,    0,    0,         cos(psi_b),      -v_b*sin(psi_b);  % s_x_dot
         0,    0,    0,         sin(psi_b),       v_b*cos(psi_b);  % s_y_dot
         0,    0,    0,         0,                0;               % delta_dot
         0,    0,    0,         0,                0;               % v_dot
         0,    0,    v_b/(l*cos_d2),  tan(delta_b)/l,  0 ];       % psi_dot

  % B_c = df/du  evaluated at (x_bar, u_bar)
  %         v_delta   a_x
  Bc = [ 0,        0;   % s_x_dot
         0,        0;   % s_y_dot
         1,        0;   % delta_dot
         0,        1;   % v_dot
         0,        0 ]; % psi_dot

  % ── Euler discretization ──────────────────────────────────────────────
  A_d = eye(5) + dt * Ac;
  B_d = dt * Bc;

  % Affine offset:  c_d = dt*(f_bar − Ac*x_bar − Bc*u_bar)
  % This guarantees  A_d*x_bar + B_d*u_bar + c_d  =  x_bar + dt*f(x_bar,u_bar)
  f_bar = ks_continuous(x_bar, u_bar, l);
  c_d   = dt * (f_bar - Ac * x_bar - Bc * u_bar);
end


% ── Local helper: continuous KS dynamics (no constraint enforcement) ──────
function f = ks_continuous(x, u, l)
  f = [ x(4)*cos(x(5));
        x(4)*sin(x(5));
        u(1);
        u(2);
        x(4)/l * tan(x(3)) ];
end
