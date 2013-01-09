%% run_optimization
% Obtain an optimized structure.

%% Description

function [z, p, state] = run_optimization(opt_prob, g, p0, options, varargin)

    p = p0;
    z = g.m(p);
    state = [];

    if ~isempty(varargin)
        state = varargin{1};
    end

    savefile = [tempdir, 'OB1__', strrep(datestr(clock), ' ', '@'), '.mat'];
    fprintf('Saving state at: %s\n', savefile);
    save_opt_state = @() save(savefile, 'k', 'opt_prob', 'g', 'z', 'p', ...
                                        'state', 'options');
    
    k = 0;
    save_opt_state();
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
        save_opt_state();
    end
end % End run_optimization function.


    
    


