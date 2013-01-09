%% run_optimization
% Obtain an optimized structure.

%% Description

function [z, p, state] = run_optimization(opt_prob, g, p0, options, varargin)

    % Set up logging.
    fname = strrep(datestr(clock), ' ', '@');
    if isempty(options.state_file)
        options.state_file = [tempdir, 'state_', fname, '.mat'];
        fprintf('Saving state data at: %s\n', options.state_file);
    end
    if isempty(options.history_file)
        options.history_file = [tempdir, 'hist_', fname, '.h5'];
        fprintf('Saving history data at: %s\n', options.history_file);
    end

    log_state = @() save(options.state_file, ...
                        'k', 'opt_prob', 'g', 'z', 'p', 'state', 'options');

    log_history = @(k, x, z, p) my_log_history(options.history_file, ...
                                        {opt_prob.get_epsilon}, k, x, z, p);
    
    p = p0;
    z = g.m(p);
    state = [];

    if ~isempty(varargin)
        state = varargin{1};
    end


    k = 0;
    log_state();

    for k = 1 : options.num_iters
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

        options.vis_progress(state, z, p);
        log_state();
        log_history(k, state.x, z, p);
    end
end % End run_optimization function.



function my_log_history(filename, get_epsilon, index, x, z, p)

    N = length(x); % Number of modes.
    addnum = @(str, num) [str, num2str(num)];

    for k = 1 : N
        get_eps{k} = @(z) my_get_eps(get_epsilon{k}, z);
        eps = get_eps{k}(z);
        dims = size(eps);
        unvec{k} = @(u) reshape(u, dims);
    end

    if index == 1 % Assume file and datasets does not yet exist, so create it.

        % Raw vector data.
        my_create_dataset(filename, 'z', length(z));
        my_create_dataset(filename, 'p', length(p));

        % Reshaped field data.
        for k = 1 : N
            my_create_dataset(filename, addnum('x', k), length(x{k}));
            my_create_dataset(filename, addnum('E', k), dims);
            my_create_dataset(filename, addnum('eps', k), dims);
        end
    end

    % Write the datasets.
    my_write_dataset(filename, 'z', z, index);
    my_write_dataset(filename, 'p', p, index);

    for k = 1 : N
        my_write_dataset(filename, addnum('x', k), x{k}, index);
        my_write_dataset(filename, addnum('E', k), unvec{k}(x{k}), index);
        my_write_dataset(filename, addnum('eps', k), get_eps{k}(z), index);
    end
end

function my_create_dataset(filename, datasetname, dims)
% Creates datasets for complex arrays.
% Includes real, imaginary, and absolute value datasets.

    type = {'_real', '_imag', '_abs'};
    for k = 1 : length(type)
        h5create(filename, ['/', datasetname, type{k}], [dims Inf], ...
                    'ChunkSize', [dims 1]);
    end
end

function my_write_dataset(filename, datasetname, data, index)
% Writes to datasets for complex arrays.
% Includes real, imaginary, and absolute value datasets.

    dims = size(data);
    if length(dims) == 2 % Detect one-dimensional arrays.
        dims = dims(1);
    end

    type = {'_real', '_imag', '_abs'};
    f = {@real, @imag, @abs};
    for k = 1 : length(type)
        h5write(filename, ['/', datasetname, type{k}], f{k}(data), ...
                [ones(1, length(dims)), index], [dims, 1]);
    end
end

function [eps] = my_get_eps(get_epsilon, z)
    epsilon = get_epsilon(z);
    eps = cat(4, epsilon{:});
end
