%% ========================================================================
%  Steady-State Convergence Analysis
%  Usage:  steadystate_convergence          (uses newest _blocks file)
%          steadystate_convergence('file')  (specific block file)
%          steadystate_convergence('file', dt)  (with timestep in s)
%  ========================================================================
function steadystate_convergence(fname, dt)
    if nargin < 1 || isempty(fname)
        % Auto-detect newest _blocks.txt file in production_output
        files = dir('../production_output/*_blocks.txt');
        if isempty(files)
            files = dir('../runtime_output/*_blocks.txt');
        end
        if isempty(files)
            % Fallback to old-style file
            if exist('../runtime_output/force_blocks.txt', 'file')
                fname = '../runtime_output/force_blocks.txt';
            elseif exist('../pipeline/force_blocks.txt', 'file')
                fname = '../pipeline/force_blocks.txt';
            else
                error('No block files found. Run a simulation first.');
            end
        else
            [~, idx] = sort([files.datenum]);
            fname = fullfile(files(idx(end)).folder, files(idx(end)).name);
        end
    end

    fprintf('Reading: %s\n', fname);

    %% 1. Read block-averaged force file
    data = fileread(fname);
    blocks = split(data, 'ITEM: TIMESTEP');
    blocks = blocks(2:end);  % discard leading empty block

    n_blocks = length(blocks);
    block_drag = zeros(n_blocks, 1);
    block_step = zeros(n_blocks, 1);

    for k = 1:n_blocks
        lines = splitlines(strtrim(blocks{k}));

        % Step number is first line
        block_step(k) = str2double(lines{1});

        % Sum Fx (col 2) across all surface triangles
        fx_sum = 0;
        for i = 5:length(lines)  % skip step + 3 header lines
            parts = strsplit(strtrim(lines{i}));
            if length(parts) >= 4
                fx_sum = fx_sum + str2double(parts{2});
            end
        end
        block_drag(k) = abs(fx_sum);
    end

    %% 2. Extract timestep
    if nargin < 2 || isempty(dt)
        % Try to get dt from the filename (production files encode params)
        % Format: force_D<diam>_V<vel>_<timestamp>_blocks.txt
        % Or try to parse it from the first block step (old format: step 1000)
        [~, name] = fileparts(fname);
        if contains(name, '_DT')
            % Extract from run_time_refinement naming
            tokens = regexp(name, 'DT([\d.e-]+)', 'tokens');
            if ~isempty(tokens)
                dt = str2double(tokens{1}{1});
            end
        end

        % If still not determined, infer from first block step
        % Standard setup: blocks every 1000 steps, so dt = step_first / 1000
        if isempty(dt) || dt <= 0
            if n_blocks >= 2
                steps_per_block = block_step(2) - block_step(1);
                % Typical dt for worst-case: ~4e-9 to 2e-8
                dt = 1e-9;  % safe default
                fprintf('  Using default dt = %.1e s (infer from block_step=%d)\n', ...
                    dt, steps_per_block);
            else
                dt = 1e-9;
            end
        end
    end
    fprintf('  Timestep dt = %.2e s\n', dt);

    %% 3. Physical time at each block (end of block window)
    phys_time = block_step * dt;  % s

    %% 4. Cumulative running average
    cum_drag = zeros(n_blocks, 1);
    for k = 1:n_blocks
        cum_drag(k) = mean(block_drag(1:k));
    end

    %% 5. Print summary
    fprintf('\nBlock  |  Step  |  Time (us)  |  Block |Fx| (N)  |  Cum. avg |Fx| (N)\n');
    fprintf('---------------------------------------------------------------------\n');
    for k = 1:n_blocks
        fprintf('  %2d   | %6d  |   %8.3f  |  %12.5e  |  %12.5e\n', ...
            k, block_step(k), phys_time(k)*1e6, block_drag(k), cum_drag(k));
    end

    %% 6. Convergence assessment
    final_drag = cum_drag(end);
    convergence_band = 0.05;  % 5% band

    % Find when cumulative average entered and stayed within ±5% of final
    within_band = abs(cum_drag - final_drag) / final_drag < convergence_band;
    first_converged = find(within_band, 1, 'first');

    fprintf('\nConvergence assessment:\n');
    fprintf('  Final cumulative drag: %.4e N\n', final_drag);
    if ~isempty(first_converged) && first_converged < n_blocks
        fprintf('  Entered ±%.0f%% band at block %d (t ≈ %.2f us)\n', ...
            convergence_band*100, first_converged, phys_time(first_converged)*1e6);
        fprintf('  ✓ Steady state reached\n');
    else
        fprintf('  ⚠ Never entered ±%.0f%% band — may not be fully converged\n', ...
            convergence_band*100);
    end

    %% 7. Plot
    figure('Color', [1 1 1], 'Position', [100 100 900, 550]);

    % Subplot 1: Block-averaged drag
    subplot(2,1,1);
    plot(phys_time*1e6, block_drag, 'o-', 'LineWidth', 1.5, ...
         'MarkerEdgeColor', [0.8500 0.3250 0.0980], ...
         'MarkerFaceColor', [0.9290 0.6940 0.1250], 'MarkerSize', 6);
    xlabel('Physical Time (\mus)');
    ylabel('Block-Averaged |F_x| (N)');
    title('Block-Averaged Drag During Steady Phase');
    grid on; grid minor;

    % Subplot 2: Cumulative running average with convergence band
    subplot(2,1,2);
    plot(phys_time*1e6, cum_drag, 's-', 'LineWidth', 2, ...
         'MarkerEdgeColor', [0 0.4470 0.7410], ...
         'MarkerFaceColor', [0.3010 0.7450 0.9330], 'MarkerSize', 6);
    xlabel('Physical Time (\mus)');
    ylabel('Cumulative Avg |F_x| (N)');
    title('Running Cumulative Average — Convergence to Steady State');
    grid on; grid minor;

    hold on;
    yline(final_drag*(1+convergence_band), '--r', sprintf('+%.0f%%', convergence_band*100), 'LineWidth', 1);
    yline(final_drag*(1-convergence_band), '--r', sprintf('-%.0f%%', convergence_band*100), 'LineWidth', 1);
    yline(final_drag, '-k', 'Final', 'LineWidth', 1);
    hold off;
end
