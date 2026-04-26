# AEROSP 740 Final Project — Linear vs. Nonlinear MPC for Vehicle Path Following

**Author:** Harsh Nikam ([nharsh@umich.edu](mailto:nharsh@umich.edu))  
**Course:** AEROSP 740 — Model Predictive Control, University of Michigan  
**Date:** April 2026

## Overview

This project compares **Linear and Nonlinear Model Predictive Control** for vehicle path following during a double lane change (DLC) maneuver. Four MPC controllers are implemented and benchmarked against a high-fidelity 29-state multi-body (MB) plant model from the [CommonRoad](https://commonroad.in.tum.de/) vehicle model library.

| Controller | Predictor Model | States | Solver |
|---|---|---|---|
| NMPC-KS | Kinematic Single-Track (KS) | 5 | CasADi / IPOPT |
| LMPC-KS | Linearized KS | 5 | MATLAB `quadprog` |
| NMPC-ST | Dynamic Single-Track (ST) + RK4 | 7 | CasADi / IPOPT |
| LMPC-ST | Linearized ST + RK4 | 7 | MATLAB `quadprog` |

The KS-based controllers replicate the reference paper by Diklic and Novoselnik (MIPRO 2024). The ST-based controllers are original contributions that introduce a dynamic bicycle model with lateral tire forces, RK4 integration, and curvature-based steady-state references.

## Key Contributions Beyond the Reference Paper

1. **Dynamic Single-Track predictor** — upgrades from 5-state kinematic to 7-state dynamic model with cornering stiffness, yaw dynamics, and side-slip.
2. **RK4 integration** — replaces forward Euler with 4th-order Runge–Kutta for more accurate state prediction inside the NLP.
3. **Curvature-based references** — physically motivated non-zero references for yaw rate and side-slip, resolving a cost-function conflict that degraded high-speed performance.

## Repository Structure

### Main Driver

| File | Description |
|---|---|
| `main_dlc_mpc_st.m` | Top-level simulation script. Runs all four controllers across multiple velocities and prediction horizons, then generates comparison plots. **Start here.** |

### MPC Controllers

| File | Description |
|---|---|
| `run_nonlinear_mpc.m` | Nonlinear MPC with KS predictor (Euler integration, CasADi/IPOPT) |
| `run_linear_mpc.m` | Linear MPC with KS predictor (analytical Jacobian, condensed QP via `quadprog`) |
| `run_nonlinear_mpc_st.m` | Nonlinear MPC with ST predictor (RK4 integration, CasADi/IPOPT) |
| `run_linear_mpc_st.m` | Linear MPC with ST predictor (numerical Jacobian, condensed QP via `quadprog`) |

### Vehicle Models and Linearization

| File | Description |
|---|---|
| `st_dynamics_continuous.m` | Continuous-time ST dynamics function (used by `linearize_st.m`) |
| `linearize_ks.m` | Analytical Jacobian of the KS model (Euler-discretized) |
| `linearize_st.m` | Numerical Jacobian of the ST model via central finite differences of the RK4 step |
| `build_prediction_matrices.m` | Builds condensed QP prediction matrices (Phi, Gamma, D) for LMPC |

### Reference and Helper Functions

| File | Description |
|---|---|
| `get_dlc_reference.m` | Computes the DLC lateral displacement and heading reference from tanh-based equations |
| `get_curvature_references.m` | Computes curvature-based steady-state references for yaw rate and side-slip |
| `simulate_plant_mb.m` | Wrapper to simulate one time step of the 29-state MB plant (via `ode45`) |
| `mb_to_ks.m` | Extracts the 5-state KS vector from the 29-state MB state |
| `mb_to_st.m` | Extracts the 7-state ST vector from the 29-state MB state |

### Plotting

| File | Description |
|---|---|
| `plot_results_st.m` | Generates all comparison figures: path following, lateral error, control inputs, and summary bar charts |
| `traj.m` | Standalone script to visualize the DLC reference path |

### External Dependency

The `commonroad-vehicle-models-master-MATLAB/` folder contains the CommonRoad vehicle model library (BMW 320i, `parameters_vehicle2.m`), which provides the 29-state multi-body plant dynamics, tire models (Pacejka PAC2002), and vehicle parameters. This library is from the [CommonRoad Vehicle Models](https://gitlab.lrz.de/tum-cps/commonroad-vehicle-models) project.

## Requirements

- **MATLAB** (R2022b or later recommended)
- **CasADi for MATLAB** — required for the NMPC controllers ([download](https://web.casadi.org/get/))
- **Optimization Toolbox** — required for `quadprog` (LMPC controllers)
- **CommonRoad vehicle models** — included in the `commonroad-vehicle-models-master-MATLAB/` folder

## How to Run

1. Open MATLAB and navigate to the project root directory.
2. Add the CommonRoad library to the MATLAB path:
   ```matlab
   addpath(genpath('commonroad-vehicle-models-master-MATLAB'));
   ```
3. Ensure CasADi is on the MATLAB path (e.g., `addpath('/path/to/casadi')`).
4. Run the main driver:
   ```matlab
   main_dlc_mpc_st
   ```
   This will execute all four controllers at velocities of 5, 10, 15, and 17 m/s with prediction horizons of 7 and 10, then produce comparison plots.

## Simulation Parameters

| Parameter | Value |
|---|---|
| Sampling time | 0.025 s |
| Simulation distance | 120 m |
| Test velocities | 5, 10, 15, 17 m/s |
| Prediction horizons | 7, 10 |
| Control horizon | Equal to prediction horizon |
| Vehicle | BMW 320i (CommonRoad vehicle2, m = 1225.887 kg) |

## Reference

V. Diklic and B. Novoselnik, "Comparison of Linear and Nonlinear Model Predictive Control for Vehicle Path Following," in *Proc. 47th MIPRO ICT and Electronics Convention*, Opatija, Croatia, May 2024, pp. 1794–1799.
