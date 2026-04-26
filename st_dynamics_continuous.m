function f = st_dynamics_continuous(x, u, p)
% ST_DYNAMICS_CONTINUOUS  Continuous-time ST (single-track) bicycle dynamics.
%
%   f = st_dynamics_continuous(x, u, p)
%
% Implements the same equations as CommonRoad vehicleDynamics_ST but
% WITHOUT constraint clamping (bounds handled by the MPC optimizer) and
% WITHOUT the low-speed kinematic fallback (handled externally if needed).
%
% State:  x = [s_x; s_y; delta; v_x; psi; psi_dot; beta]   (7 x 1)
% Input:  u = [v_delta; a_x]                                 (2 x 1)
%
% The dynamics include:
%   - Lateral tire forces via linear cornering stiffness from Pacejka p_ky1/p_dy1
%   - Load transfer due to longitudinal acceleration
%   - Yaw moment balance  (I_z * psi_ddot = ...)
%   - Side-slip rate      (beta_dot = ...)
%
% Inputs
%   x  - State vector (7 x 1)
%   u  - Input vector (2 x 1)
%   p  - Vehicle parameter struct from parameters_vehicle2()
%
% Output
%   f  - Time derivative dx/dt (7 x 1)

  % Physical constants
  g = 9.81;  % [m/s^2]

  % Vehicle parameters
  mu  = p.tire.p_dy1;                     % friction coefficient
  C_Sf = -p.tire.p_ky1 / p.tire.p_dy1;   % front cornering stiffness
  C_Sr = -p.tire.p_ky1 / p.tire.p_dy1;   % rear cornering stiffness
  lf  = p.a;                              % CG to front axle [m]
  lr  = p.b;                              % CG to rear axle  [m]
  lwb = lf + lr;                          % wheelbase [m]
  h   = p.h_s;                            % sprung mass CG height [m]
  m   = p.m;                              % vehicle mass [kg]
  I   = p.I_z;                            % yaw inertia [kg m^2]

  % Unpack state
  % s_x   = x(1);  % not used in dynamics
  % s_y   = x(2);  % not used in dynamics
  delta   = x(3);
  v_x     = x(4);
  psi     = x(5);
  psi_dot = x(6);
  beta    = x(7);

  % Unpack input
  v_delta = u(1);
  a_x     = u(2);

  % Guard against very low speed (avoid division by zero)
  % At low speed, fall back to kinematic model
  if abs(v_x) < 0.5
    % Kinematic model (COG reference) for very low speed
    beta_kin = atan(tan(delta) * lr / lwb);
    f = zeros(7, 1);
    f(1) = v_x * cos(beta_kin + psi);
    f(2) = v_x * sin(beta_kin + psi);
    f(3) = v_delta;
    f(4) = a_x;
    f(5) = v_x * cos(beta_kin) * tan(delta) / lwb;
    f(6) = 0;   % psi_ddot ≈ 0 at near-zero speed
    f(7) = 0;   % beta_dot ≈ 0 at near-zero speed
    return;
  end

  % Gravity + load-transfer terms (NOT actual forces — matches CommonRoad
  % vehicleDynamics_ST which uses these raw terms directly in the equations).
  %   glr = g*lr - a_x*h   (front-axle related term)
  %   glf = g*lf + a_x*h   (rear-axle related term)
  glr = g * lr - a_x * h;    % [m^2/s^2]  (front)
  glf = g * lf + a_x * h;    % [m^2/s^2]  (rear)

  % Dynamic single-track equations — identical to CommonRoad vehicleDynamics_ST
  f = zeros(7, 1);
  f(1) = v_x * cos(beta + psi);                                % dot_s_x
  f(2) = v_x * sin(beta + psi);                                % dot_s_y
  f(3) = v_delta;                                               % dot_delta
  f(4) = a_x;                                                  % dot_v_x
  f(5) = psi_dot;                                               % dot_psi
  f(6) = -mu*m/(v_x*I*lwb) * ...                                % dot_psi_dot
         (lf^2*C_Sf*glr + lr^2*C_Sr*glf) * psi_dot ...
       + mu*m/(I*lwb) * ...
         (lr*C_Sr*glf - lf*C_Sf*glr) * beta ...
       + mu*m/(I*lwb) * lf*C_Sf*glr * delta;
  f(7) = (mu/(v_x^2*lwb) * ...                                  % dot_beta
         (C_Sr*glf*lr - C_Sf*glr*lf) - 1) * psi_dot ...
       - mu/(v_x*lwb) * ...
         (C_Sr*glf + C_Sf*glr) * beta ...
       + mu/(v_x*lwb) * C_Sf*glr * delta;
end
