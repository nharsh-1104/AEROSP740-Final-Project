%% main_dlc_mpc_st.m
% =========================================================================
%  Extended Comparison: KS-based vs ST-based MPC for Vehicle Path Following
% =========================================================================
%  Plant:       CommonRoad Multi-Body (MB) Model 2 (BMW 320i, 29 states)
%
%  Predictors:
%    KS — Kinematic Single-Track  (5 states, Euler discretization)
%    ST — Dynamic Single-Track    (7 states, RK4 discretization)
%         Includes cornering stiffness, yaw rate, side-slip dynamics
%
%  Controllers:
%    1. NMPC_KS  — Nonlinear MPC  with KS predictor (original paper)
%    2. LMPC_KS  — Linear MPC     with KS predictor (original paper)
%    3. NMPC_ST  — Nonlinear MPC  with ST predictor + RK4  (improved)
%    4. LMPC_ST  — Linear MPC     with ST predictor + RK4  (improved)
%
%  Maneuver:   Double Lane Change (DLC) per Diklic & Novoselnik, MIPRO 2024
%
%  Dependencies: CasADi, MATLAB Optimization Toolbox
% =========================================================================

clear; clc; close all;

%% ── 1.  Paths ──────────────────────────────────────────────────────────
base    = fileparts(mfilename('fullpath'));
cr_root = fullfile(base, 'commonroad-vehicle-models-master-MATLAB', 'MATLAB');

addpath(fullfile(cr_root, 'vehiclemodels'));
addpath(fullfile(cr_root, 'vehiclemodels', 'utils'));
addpath(fullfile(cr_root, 'vehiclemodels', 'utils', 'tireModel'));
addpath(fullfile(cr_root, 'vehiclemodels', 'utils', 'unitConversions'));

import casadi.*;

%% ── 2.  Vehicle parameters ────────────────────────────────────────────
p = parameters_vehicle2();
l = 2.5;    % wheelbase for KS bicycle model [m] (paper value)

fprintf('Vehicle: BMW 320i (CommonRoad vehicle2)\n');
fprintf('  Wheelbase (KS)  l = %.3f m\n', l);
fprintf('  Wheelbase (ST)  l = %.3f m  (a=%.3f + b=%.3f)\n', p.a+p.b, p.a, p.b);
fprintf('  Mass  m = %.1f kg,  I_z = %.1f kg m^2\n', p.m, p.I_z);

%% ── 3.  DLC parameters (Table I of the paper) ─────────────────────────
dlc.shape = 2.4;
dlc.dx1   = 25.00;    dlc.dx2 = 21.95;
dlc.dy1   =  4.05;    dlc.dy2 =  5.70;
dlc.Xs1   = 27.19;    dlc.Xs2 = 56.46;

%% ── 4.  MPC settings ──────────────────────────────────────────────────
mpc.dt    = 0.025;    % sample time [s] (paper value)
mpc.nu    = 2;

% ── KS (5-state) settings ──────────────────────────────────────────────
mpc.nx    = 5;
mpc.l     = l;
%   Paper weights reordered for KS state [s_x, s_y, delta, v, psi]
mpc.Q     = diag([ 20,   20,  0.1,   5,  200 ]);
mpc.R     = diag([ 1,    1 ]);
mpc.QN    = diag([ 100, 100,  0.1,   5, 1000 ]);

% ── ST (7-state) settings ──────────────────────────────────────────────
mpc.nx_st = 7;
%   States: [s_x, s_y, delta, v_x, psi, psi_dot, beta]
%   Same position/heading/velocity weights as KS, plus weights on the two
%   new states: psi_dot and beta.  Penalising these encourages stability.
mpc.Q_st  = diag([ 20,   20,  0.1,   5,  200,   10,   10 ]);
mpc.QN_st = diag([ 100, 100,  0.1,   5, 1000,   50,   50 ]);

% Input bounds (same for both)
mpc.u_min = [p.steering.v_min; -p.longitudinal.a_max];
mpc.u_max = [p.steering.v_max;  p.longitudinal.a_max];
mpc.d_min = p.steering.min;
mpc.d_max = p.steering.max;

%% ── 5.  Simulation sweep ──────────────────────────────────────────────
velocities = [5, 10, 15, 17];
hp_values  = [7, 10];        % focus on medium/long horizons for ST benefit
sim_dist   = 120;             % [m]

results   = struct();
sim_count = 0;

for vi = 1:numel(velocities)
  v0    = velocities(vi);
  N_sim = ceil(sim_dist / (v0 * mpc.dt));

  for hi = 1:numel(hp_values)
    hp = hp_values(hi);
    hu = hp;   % use full control horizon for cleaner comparison

    sim_count = sim_count + 1;

    fprintf('\n%s\n  v0 = %d m/s | hp = %d | hu = %d\n%s\n', ...
            repmat('=',1,50), v0, hp, hu, repmat('=',1,50));

    % ── Initial conditions ───────────────────────────────────────────
    x_ks0 = [0; 0; 0; v0; 0];
    x_st0 = [0; 0; 0; v0; 0; 0; 0];
    x_mb0 = init_MB([0; 0; 0; v0; 0; 0; 0], p)';

    % ── 1. NMPC with KS predictor ───────────────────────────────────
    fprintf('  NMPC_KS ... ');
    try
      r1 = run_nonlinear_mpc(x_ks0, x_mb0, N_sim, hp, hu, mpc, dlc, p, v0);
      fprintf('eps=%.5f  tc=%.4fs\n', r1.eps, r1.t_avg);
    catch ME, fprintf('FAILED: %s\n', ME.message); r1=[]; end

    % ── 2. LMPC with KS predictor ───────────────────────────────────
    fprintf('  LMPC_KS ... ');
    try
      r2 = run_linear_mpc(x_ks0, x_mb0, N_sim, hp, hu, mpc, dlc, p, v0);
      fprintf('eps=%.5f  tc=%.4fs\n', r2.eps, r2.t_avg);
    catch ME, fprintf('FAILED: %s\n', ME.message); r2=[]; end

    % ── 3. NMPC with ST predictor + RK4 ─────────────────────────────
    fprintf('  NMPC_ST ... ');
    try
      r3 = run_nonlinear_mpc_st(x_st0, x_mb0, N_sim, hp, hu, mpc, dlc, p, v0);
      fprintf('eps=%.5f  tc=%.4fs\n', r3.eps, r3.t_avg);
    catch ME, fprintf('FAILED: %s\n', ME.message); r3=[]; end

    % ── 4. LMPC with ST predictor ───────────────────────────────────
    fprintf('  LMPC_ST ... ');
    try
      r4 = run_linear_mpc_st(x_st0, x_mb0, N_sim, hp, hu, mpc, dlc, p, v0);
      fprintf('eps=%.5f  tc=%.4fs\n', r4.eps, r4.t_avg);
    catch ME, fprintf('FAILED: %s\n', ME.message); r4=[]; end

    results(sim_count).v0      = v0;
    results(sim_count).hp      = hp;
    results(sim_count).hu      = hu;
    results(sim_count).nmpc_ks = r1;
    results(sim_count).lmpc_ks = r2;
    results(sim_count).nmpc_st = r3;
    results(sim_count).lmpc_st = r4;
  end
end

%% ── 6.  Summary table ─────────────────────────────────────────────────
fprintf('\n\n');
fprintf('%-5s %-4s | %-12s %-10s | %-12s %-10s | %-12s %-10s | %-12s %-10s\n', ...
  'v0','hp', 'NMPC_KS_e','tc_KS_NL', 'LMPC_KS_e','tc_KS_L', ...
              'NMPC_ST_e','tc_ST_NL', 'LMPC_ST_e','tc_ST_L');
fprintf('%s\n', repmat('-',1,115));

for i = 1:sim_count
  r = results(i);
  vals = {'N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A'};
  flds = {'nmpc_ks','lmpc_ks','nmpc_st','lmpc_st'};
  for fi = 1:4
    if ~isempty(r.(flds{fi}))
      vals{2*fi-1} = sprintf('%.5f', r.(flds{fi}).eps);
      vals{2*fi}   = sprintf('%.4f',  r.(flds{fi}).t_avg);
    end
  end
  fprintf('%-5d %-4d | %-12s %-10s | %-12s %-10s | %-12s %-10s | %-12s %-10s\n', ...
          r.v0, r.hp, vals{:});
end

save('mpc_results_st.mat', 'results');
fprintf('\nResults saved to mpc_results_st.mat\n');

%% ── 7.  Plots ─────────────────────────────────────────────────────────
plot_results_st(results, dlc);
