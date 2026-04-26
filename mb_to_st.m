function x_st = mb_to_st(x_mb)
% MB_TO_ST  Extract the 7 ST-observable states from a 29-state MB vector.
%
%   x_st = mb_to_st(x_mb)
%
% The ST (single-track) predictor model uses state:
%   x = [s_x, s_y, delta, v_x, psi, psi_dot, beta]
%
% Mapping from MB state vector:
%   MB x(1)  = s_x        -> ST x(1)
%   MB x(2)  = s_y        -> ST x(2)
%   MB x(3)  = delta      -> ST x(3)
%   MB x(4)  = v_x        -> ST x(4)   (longitudinal velocity)
%   MB x(5)  = psi        -> ST x(5)   (yaw angle)
%   MB x(6)  = psi_dot    -> ST x(6)   (yaw rate)
%   MB x(11) = v_y        -> used to compute beta
%
%   beta = atan2(v_y, v_x)             -> ST x(7)   (side-slip angle)
%
% Input
%   x_mb  - MB state vector (29 x 1)
%
% Output
%   x_st  - ST state vector (7 x 1)

  v_x = x_mb(4);
  v_y = x_mb(11);

  % Side-slip angle at centre of mass
  if abs(v_x) < 0.01
    beta = 0;
  else
    beta = atan2(v_y, v_x);
  end

  x_st = [x_mb(1);   % s_x
          x_mb(2);   % s_y
          x_mb(3);   % delta
          x_mb(4);   % v_x
          x_mb(5);   % psi
          x_mb(6);   % psi_dot
          beta];     % beta
end
