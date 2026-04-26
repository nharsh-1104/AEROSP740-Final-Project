function x_next = simulate_plant_mb(x_mb, u, p, dt)
% SIMULATE_PLANT_MB  Integrate the CommonRoad multi-body plant one step.
%
%   x_next = simulate_plant_mb(x_mb, u, p, dt)
%
% Uses ode45 to integrate vehicleDynamics_MB from t = 0 to t = dt.
%
% Inputs
%   x_mb  - Current MB state vector  (29 x 1)
%            Layout (see init_MB.m / vehicleDynamics_MB.m):
%              x(1)  = s_x        x(2)  = s_y
%              x(3)  = delta      x(4)  = v_x
%              x(5)  = psi        x(6)  = psi_dot
%              x(7)  = phi_s      x(8)  = phi_s_dot
%              x(9)  = theta_s    x(10) = theta_s_dot
%              x(11) = v_y        x(12) = z_s
%              x(13) = w_s        (... 16 more suspension/wheel states)
%   u     - Control input [v_delta (rad/s); a_x (m/s^2)]  (2 x 1)
%   p     - Vehicle parameter struct from parameters_vehicle2()
%   dt    - Integration step [s]
%
% Output
%   x_next - MB state at t = dt  (29 x 1)
%
% Notes
%  The vehicleDynamics_MB function already enforces steering and
%  acceleration constraints internally via steeringConstraints and
%  accelerationConstraints.  No additional saturation is applied here.

  odefun = @(~, x) vehicleDynamics_MB(x, u, p);

  opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-6 * ones(29,1));  % paper: rtol=1e-3, atol=1e-6
  [~, X] = ode45(odefun, [0, dt], x_mb(:), opts);
  x_next = X(end, :)';
end
