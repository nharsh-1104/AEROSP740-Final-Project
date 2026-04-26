function [Phi, Gamma, D] = build_prediction_matrices(A_d, B_d, c_d, nx, nu, hp, hu)
% BUILD_PREDICTION_MATRICES  Condensed prediction matrices for linear MPC.
%
%   [Phi, Gamma, D] = build_prediction_matrices(A_d, B_d, c_d, nx, nu, hp, hu)
%
% Given the discrete affine model  x_{k+1} = A_d*x_k + B_d*u_k + c_d,
% this function constructs matrices such that:
%
%   X_pred  =  Phi * x_0  +  Gamma * U  +  D
%
% where
%   X_pred  = [x_1; x_2; ...; x_hp]           (hp*nx × 1)
%   U       = [u_0; u_1; ...; u_{hu-1}]        (hu*nu × 1)   decision vars
%
% Control-horizon convention: u_k = u_{hu-1}  for all k ≥ hu.
%
% Inputs
%   A_d  - Discrete state matrix   (nx × nx)
%   B_d  - Discrete input matrix   (nx × nu)
%   c_d  - Discrete affine offset  (nx × 1)
%   nx   - State dimension
%   nu   - Input dimension
%   hp   - Prediction horizon
%   hu   - Control horizon  (hu ≤ hp)
%
% Outputs
%   Phi    - Initial-state propagation matrix  (hp*nx × nx)
%   Gamma  - Input sensitivity matrix         (hp*nx × hu*nu)
%   D      - Bias (affine offset) vector      (hp*nx × 1)

  %% ── Phi and D (do not depend on hu) ──────────────────────────────────
  Phi = zeros(hp * nx, nx);
  D   = zeros(hp * nx, 1);

  A_pow = eye(nx);     % will hold A_d^k
  D_k   = zeros(nx, 1);  % D_0 = 0

  for k = 1:hp
    A_pow = A_d * A_pow;         % A_d^k
    D_k   = A_d * D_k + c_d;    % sum_{j=0}^{k-1} A_d^j * c_d

    rows = (k-1)*nx + (1:nx);
    Phi(rows, :) = A_pow;
    D(rows)      = D_k;
  end

  %% ── Full Gamma_full (hp*nx × hp*nu) ──────────────────────────────────
  % Gamma_full(k, j) = A_d^{k-j} * B_d,   for 1-indexed j = 1..k
  % (corresponds to u_{j-1} acting on state at prediction step k)
  Gamma_full = zeros(hp * nx, hp * nu);

  for k = 1:hp
    rows = (k-1)*nx + (1:nx);
    % Iterate from j=k down to j=1; accumulate A_d powers from 0 upward
    A_pow_j = eye(nx);      % A_d^(k-j) when j=k  →  A_d^0
    for j = k:-1:1
      cols = (j-1)*nu + (1:nu);
      Gamma_full(rows, cols) = A_pow_j * B_d;
      A_pow_j = A_d * A_pow_j;  % prepare A_d^(k-j+1) for next j
    end
  end

  %% ── Reduction matrix S: (hp*nu × hu*nu) ─────────────────────────────
  % Maps compressed decision vector U (hu inputs) to full input sequence
  % (hp inputs) by repeating the last decision variable for k ≥ hu.
  S = zeros(hp * nu, hu * nu);
  for j = 1:hp
    j_hu = min(j, hu);
    S((j-1)*nu+(1:nu), (j_hu-1)*nu+(1:nu)) = eye(nu);
  end

  %% ── Reduced Gamma ────────────────────────────────────────────────────
  Gamma = Gamma_full * S;    % (hp*nx × hu*nu)
end
