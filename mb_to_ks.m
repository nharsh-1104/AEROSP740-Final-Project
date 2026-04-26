function x_ks = mb_to_ks(x_mb)
% MB_TO_KS  Extract the 5 KS-observable states from a 29-state MB vector.
%
%   x_ks = mb_to_ks(x_mb)
%
% The KS predictor model uses state  x = [s_x, s_y, delta, v_x, psi].
% These correspond directly to the first 5 elements of the MB state vector:
%
%   MB x(1) = s_x    (global x-position)    →  KS x(1)
%   MB x(2) = s_y    (global y-position)    →  KS x(2)
%   MB x(3) = delta  (front steering angle) →  KS x(3)
%   MB x(4) = v_x    (longitudinal velocity)→  KS x(4)
%   MB x(5) = psi    (yaw angle)            →  KS x(5)
%
% The remaining 24 MB states — yaw rate, roll, pitch, lateral velocity,
% suspension deflections, wheel angular speeds, and compliance terms —
% are NOT visible to the KS predictor.  This structural model-plant
% mismatch is central to the comparison study.
%
% Input
%   x_mb  - MB state vector (29 x 1)
%
% Output
%   x_ks  - KS state vector  (5 x 1)

  x_ks = x_mb(1:5);
end
