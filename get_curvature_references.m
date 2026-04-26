function [psi_dot_ref, beta_ref] = get_curvature_references(sx, dlc, v0, lr)
% GET_CURVATURE_REFERENCES  Steady-state psi_dot and beta references from path curvature.
%
%   [psi_dot_ref, beta_ref] = get_curvature_references(sx, dlc, v0, lr)
%
% Computes the path curvature kappa = dpsi/dsx at each longitudinal
% position, then derives the steady-state yaw rate and side-slip angle
% that a bicycle model would need to track the DLC path:
%
%   psi_dot_ref = v0 * kappa
%   beta_ref    = lr * kappa      (rear-axle kinematic approximation)
%
% Inputs
%   sx   - Longitudinal positions (1 x N or N x 1)
%   dlc  - DLC parameter struct
%   v0   - Reference longitudinal speed [m/s]
%   lr   - Distance from CG to rear axle [m]  (= p.b)
%
% Outputs
%   psi_dot_ref - Reference yaw rate     [rad/s]  (same size as sx)
%   beta_ref    - Reference side-slip    [rad]     (same size as sx)

  n = numel(sx);

  % Get psi_ref from the DLC path at a finer grid for accurate derivatives
  [~, psi_ref] = get_dlc_reference(sx, dlc);

  % Curvature via central finite differences on psi_ref vs sx
  kappa = zeros(size(sx));

  if n >= 3
    % Central differences for interior points
    for i = 2:n-1
      ds = sx(i+1) - sx(i-1);
      if abs(ds) > 1e-10
        kappa(i) = (psi_ref(i+1) - psi_ref(i-1)) / ds;
      end
    end
    % Forward difference for first point
    ds1 = sx(2) - sx(1);
    if abs(ds1) > 1e-10
      kappa(1) = (psi_ref(2) - psi_ref(1)) / ds1;
    end
    % Backward difference for last point
    dsn = sx(n) - sx(n-1);
    if abs(dsn) > 1e-10
      kappa(n) = (psi_ref(n) - psi_ref(n-1)) / dsn;
    end
  elseif n == 2
    ds = sx(2) - sx(1);
    if abs(ds) > 1e-10
      kappa(:) = (psi_ref(2) - psi_ref(1)) / ds;
    end
  end

  % Steady-state references from bicycle model
  psi_dot_ref = v0 * kappa;
  beta_ref    = lr * kappa;
end
