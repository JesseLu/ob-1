%% line_search_brute
% Search for the minimum along a line function.

%% Description
% In essence, just continuous trying points along a line until
% a minimum is found.

function [optim_step] = line_search_convex(f, dx, x0, step_err)


    % Search for an increase in f.
    step = 1;
    step_growth_rate = 1.1;

    f_min = f(x0);
    step_min = 0;
    step_max = 1;
    for k = 1 : 1e3 % Set max to 1000 steps.
        f_max = f(x0 + step_max * dx);

        if f_max > f_min % We have a maximum bound for the step length.
            break
        else
            step_max = step_max * step_growth_rate;
        end
    end

    % "Ten-sect" to find the optimal step to within step_err.
    norm_dx = norm(dx);
    while (step_max - step_min) > step_err/norm_dx
        ds = (step_max - step_min) / 10;
        new_steps = [step_min:ds:step_max, step_max];
        
        for l = 1 : length(new_steps)
            new_f(l) = f(x0 + new_steps(l) * dx);
        end

        [f_temp, ind] = min(new_f);
        
        step_min = new_steps(max([ind-1, 1]));
        step_max = new_steps(min([ind+1, length(new_steps)]));
    end

    % Return optimal step.
    optim_step = new_steps(ind);

