function [y_ref, psi_ref] = get_dlc_reference(px, dlc)
% GET_DLC_REFERENCE Compute DLC reference lateral position and yaw angle
%
% Implements equations (1)-(4) from the paper:
%   y_ref   = (dy1/2)*(1+tanh(z1)) - (dy2/2)*(1+tanh(z2))
%   psi_ref = arctan( (dy1*1.2)/(dx1*cosh(z1)^2) - (dy2*1.2)/(dx2*cosh(z2)^2) )
%
% Inputs:
%   px  - Current longitudinal position (scalar or vector) [m]
%   dlc - Struct with DLC parameters:
%           .shape, .dx1, .dx2, .dy1, .dy2, .Xs1, .Xs2
%
% Outputs:
%   y_ref   - Reference lateral position [m]
%   psi_ref - Reference yaw angle [rad]

    % Compute z1 and z2 (eq. 3 and 4)
    z1 = (dlc.shape ./ dlc.dx1) .* (px - dlc.Xs1) - (dlc.shape / 2);
    z2 = (dlc.shape ./ dlc.dx2) .* (px - dlc.Xs2) - (dlc.shape / 2);

    % Reference lateral position (eq. 1)
    y_ref = (dlc.dy1 / 2) .* (1 + tanh(z1)) - (dlc.dy2 / 2) .* (1 + tanh(z2));

    % Reference yaw angle (eq. 2): psi_ref = atan(dy/dx)
    %   dy/dx = (dy1*shape)/(2*dx1*cosh(z1)^2) - (dy2*shape)/(2*dx2*cosh(z2)^2)
    %   shape/2 = 1.2 when shape = 2.4
    half_shape = dlc.shape / 2;
    psi_ref = atan( ...
        (dlc.dy1 * half_shape) ./ (dlc.dx1 .* cosh(z1).^2) - ...
        (dlc.dy2 * half_shape) ./ (dlc.dx2 .* cosh(z2).^2) );

end
