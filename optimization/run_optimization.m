%% run_optimization
% Obtain an optimized structure.

%% Description

function [z, p, state] = run_optimization(opt_prob, g, p0, options, varargin)

    %% Set up logging.
    options = check_file_names(options);
    log_state = @() save(options.state_file, ...
                        'k', 'opt_prob', 'g', 'z', 'p', ...
                        'state', 'options', 'progress');
    log_history = @(k, x, z, p) history_logger(options.history_file, ...
                                        {opt_prob.get_epsilon}, k, x, z, p);
    

    %% Initialize variables
    p = p0;
    z = g.m(p);
    state = [];
    progress = [];

    if ~isempty(varargin)
        state = varargin{1};
        progress = varargin{2};
    end

    k = options.starting_iter - 1;
    log_state(); % Log initial state.
    vp_state = [];


    %% Run optimization
    for k = options.starting_iter : options.num_iters
        fprintf('%d: ', k);

        % Generate Q(z).
        if strcmp(options.paradigm, 'local')
            [P, q, state] = prodQ_local(z, opt_prob, state, ...
                                        options.paradigm_args{:});
                                    
        elseif strcmp(options.paradigm, 'global')
            [P, q, state] = prodQ_global(z, opt_prob, state, ...
                                        options.paradigm_args{:});
        else
            error('Invalid paradigm.');
        end

        % Update the structure variable.
        [z, p] = update_structure(P, q, g, p, options.structure_args{:});  

        % Log and visualize.
        progress = options.vis_progress(k, state.x, z, p);
        log_state();
        log_history(k, state.x, z, p); % Do this last, in case hdf5 calls fail.

        fprintf('\n');
    end
end % End run_optimization function.


%% Private functions
function [options] = check_file_names(options)
% Generate default names if needed for the state and history files.
    fname = strrep(datestr(clock), ' ', '@');
    if isempty(options.state_file)
        options.state_file = [tempdir, 'state_', fname, '.mat'];
        fprintf('Saving state data at: %s\n', options.state_file);
    end
    if isempty(options.history_file)
        options.history_file = [tempdir, 'hist_', fname, '.h5'];
        fprintf('Saving history data at: %s\n', options.history_file);
    end
end

