%% run_optimization
% Obtain an optimized structure.

%% Description

function [z, p, state] = run_optimization(opt_prob, g, p0, paradigm)

    %% Parse inputs

    p = p0;
    z = g.m(p);

    % Form initial state

    state = [];

    for k = 1 : 10
        if strcmp(paradigm, 'local')
            [P, q, state] = prodQ_local(z, opt_prob, state);
        elseif strcmp(paradigm, 'global')
            [P, q, state] = prodQ_global(z, opt_prob, state);
        else
            error('Invalid paradigm.');
        end

        % Update the structure variable.
        [z, p] = update_structure(P, q, g, p);  
    end



