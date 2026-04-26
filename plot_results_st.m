function plot_results_st(results, dlc)
% PLOT_RESULTS_ST  Comparison plots: KS-based vs ST-based MPC controllers.
%
%   plot_results_st(results, dlc)
%
% Produces per-velocity figures:
%   1. Path following  (sy vs sx)  — all 4 controllers + reference
%   2. Lateral error   (sy - sy_ref vs sx)
%   3. Control inputs  (steering velocity & longitudinal acceleration)
%   4. Bar chart of error (eps) and computation time across all cases

  velocities = unique([results.v0]);
  hp_vals    = unique([results.hp]);
  n_v  = numel(velocities);
  n_hp = numel(hp_vals);

  % Controller field names and display properties
  ctrl = struct( ...
    'field', {'nmpc_ks', 'lmpc_ks', 'nmpc_st', 'lmpc_st'}, ...
    'label', {'NMPC (KS)', 'LMPC (KS)', 'NMPC (ST+RK4)', 'LMPC (ST+RK4)'}, ...
    'color', {[0.00 0.45 0.74], [0.85 0.33 0.10], ...     % blue, orange
              [0.47 0.67 0.19], [0.64 0.08 0.18]}, ...     % green, crimson
    'style', {'-', '-', '--', '--'} );

  %% ── Figure set 1: Path following ─────────────────────────────────────
  for vi = 1:n_v
    v0 = velocities(vi);
    figure('Name', sprintf('Path — v0=%d m/s', v0), ...
           'NumberTitle','off', 'Position', [80 80 950 300*n_hp]);
    sgtitle(sprintf('DLC Path Following  (v_0 = %d m/s,  MB plant)', v0), ...
            'FontSize', 13, 'FontWeight', 'bold');

    for hi = 1:n_hp
      hp  = hp_vals(hi);
      idx = find([results.v0]==v0 & [results.hp]==hp, 1);
      if isempty(idx), continue; end
      r = results(idx);

      subplot(n_hp,1,hi); hold on; grid on; box on;

      % Reference
      sx_any = get_first_sx(r, ctrl);
      if isempty(sx_any), continue; end
      [sy_ref_p, ~] = get_dlc_reference(sx_any, dlc);
      plot(sx_any, sy_ref_p, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference');

      % Each controller
      for ci = 1:4
        ri = r.(ctrl(ci).field);
        if isempty(ri), continue; end
        plot(ri.sx, ri.sy, ctrl(ci).style, 'Color', ctrl(ci).color, ...
             'LineWidth', 1.8, ...
             'DisplayName', sprintf('%s (\\epsilon=%.4f)', ctrl(ci).label, ri.eps));
      end

      xlabel('s_x [m]'); ylabel('s_y [m]');
      title(sprintf('h_p = %d', hp));
      legend('Location','best','FontSize',7);
      xlim([0 max(sx_any)]);
    end
  end

  %% ── Figure set 2: Lateral error ──────────────────────────────────────
  for vi = 1:n_v
    v0 = velocities(vi);
    figure('Name', sprintf('Lat.Error — v0=%d m/s', v0), ...
           'NumberTitle','off', 'Position', [120 120 950 300*n_hp]);
    sgtitle(sprintf('Lateral Tracking Error  (v_0 = %d m/s)', v0), ...
            'FontSize', 13, 'FontWeight', 'bold');

    for hi = 1:n_hp
      hp  = hp_vals(hi);
      idx = find([results.v0]==v0 & [results.hp]==hp, 1);
      if isempty(idx), continue; end
      r = results(idx);

      subplot(n_hp,1,hi); hold on; grid on; box on;

      for ci = 1:4
        ri = r.(ctrl(ci).field);
        if isempty(ri), continue; end
        err = ri.sy - ri.sy_ref;
        plot(ri.sx, err, ctrl(ci).style, 'Color', ctrl(ci).color, ...
             'LineWidth', 1.5, ...
             'DisplayName', sprintf('%s (\\epsilon=%.4f)', ctrl(ci).label, ri.eps));
      end
      yline(0,'k--','LineWidth',0.8);
      xlabel('s_x [m]'); ylabel('Lateral error [m]');
      title(sprintf('h_p = %d', hp));
      legend('Location','best','FontSize',7);
    end
  end

  %% ── Figure set 3: Control inputs (steering rate + acceleration) ──────
  for vi = 1:n_v
    v0 = velocities(vi);
    % Use first available hp for control plots
    hp = hp_vals(1);
    idx = find([results.v0]==v0 & [results.hp]==hp, 1);
    if isempty(idx), continue; end
    r = results(idx);

    figure('Name', sprintf('Controls — v0=%d', v0), ...
           'NumberTitle','off', 'Position', [160 160 950 500]);
    sgtitle(sprintf('Control Inputs  (v_0 = %d m/s,  h_p = %d)', v0, hp), ...
            'FontSize', 13, 'FontWeight', 'bold');

    % Steering velocity
    subplot(2,1,1); hold on; grid on; box on;
    for ci = 1:4
      ri = r.(ctrl(ci).field);
      if isempty(ri), continue; end
      steps = 0:size(ri.u,2)-1;
      plot(steps, ri.u(1,:)*180/pi, ctrl(ci).style, ...
           'Color', ctrl(ci).color, 'LineWidth', 1.5, ...
           'DisplayName', ctrl(ci).label);
    end
    xlabel('step'); ylabel('v_{\delta} [deg/s]');
    title('Steering velocity'); legend('Location','best','FontSize',7);

    % Longitudinal acceleration
    subplot(2,1,2); hold on; grid on; box on;
    for ci = 1:4
      ri = r.(ctrl(ci).field);
      if isempty(ri), continue; end
      plot(ri.u(2,:), ctrl(ci).style, 'Color', ctrl(ci).color, ...
           'LineWidth', 1.5, 'DisplayName', ctrl(ci).label);
    end
    yline(0,'k--','LineWidth',0.8);
    xlabel('step'); ylabel('a_x [m/s^2]');
    title('Longitudinal acceleration'); legend('Location','best','FontSize',7);
  end

  %% ── Figure 4: Bar chart — error & computation time ───────────────────
  figure('Name','Summary','NumberTitle','off','Position',[200 200 1100 500]);

  labels = {};
  eps_data = zeros(0, 4);
  tc_data  = zeros(0, 4);

  for vi = 1:n_v
    v0 = velocities(vi);
    for hi = 1:n_hp
      hp = hp_vals(hi);
      idx = find([results.v0]==v0 & [results.hp]==hp, 1);
      if isempty(idx), continue; end
      r = results(idx);
      labels{end+1} = sprintf('v=%d, h_p=%d', v0, hp);  %#ok<AGROW>
      row_e = NaN(1,4);  row_t = NaN(1,4);
      for ci = 1:4
        ri = r.(ctrl(ci).field);
        if ~isempty(ri)
          row_e(ci) = ri.eps;
          row_t(ci) = ri.t_avg;
        end
      end
      eps_data(end+1,:) = row_e;  %#ok<AGROW>
      tc_data(end+1,:)  = row_t;  %#ok<AGROW>
    end
  end

  x_pos = 1:numel(labels);

  subplot(1,2,1);
  bh = bar(x_pos, eps_data, 'grouped');
  for ci = 1:4, bh(ci).FaceColor = ctrl(ci).color; end
  set(gca,'XTick',x_pos,'XTickLabel',labels,'XTickLabelRotation',30);
  ylabel('\epsilon  (position MSE)'); title('Tracking Error');
  legend({ctrl.label},'Location','northwest','FontSize',7);
  grid on; box on;

  subplot(1,2,2);
  bh2 = bar(x_pos, tc_data, 'grouped');
  for ci = 1:4, bh2(ci).FaceColor = ctrl(ci).color; end
  set(gca,'XTick',x_pos,'XTickLabel',labels,'XTickLabelRotation',30);
  ylabel('Average solve time [s]'); title('Computation Time');
  legend({ctrl.label},'Location','northwest','FontSize',7);
  grid on; box on;
end


% ═══════════════════════════════════════════════════════════════════════════
%  Local helpers
% ═══════════════════════════════════════════════════════════════════════════
function sx = get_first_sx(r, ctrl)
% Return the sx vector from the first non-empty controller result.
  sx = [];
  for ci = 1:numel(ctrl)
    ri = r.(ctrl(ci).field);
    if ~isempty(ri), sx = ri.sx; return; end
  end
end
