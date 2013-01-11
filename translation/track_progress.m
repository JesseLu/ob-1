%% Private functions
function [progress_out, x] = track_progress(opt_prob, struct_obj, mode_sel, ...
                                        k, x, z, p, varargin)

    N = length(opt_prob);

    persistent progress
    if (k == 1)  
        for i = 1 : N
            empty_cells{i} = [];
        end
        progress = struct(  'iters', [], ... 
                            'out_power', {empty_cells}, ...
                            'out_degrees', {empty_cells}, ...
                            'res_norm', [], ...
                            'struct_obj', []);
    end

    % Run verification layer to calculate relevant physical quantities
    if isempty(x)
        modes = verification_layer(opt_prob, z);
    else
        modes = verification_layer(opt_prob, z, x);
    end

    progress.iters(end+1) = k;
    progress.struct_obj(k) = struct_obj.w(p);
    for i = 1 : N
        progress.out_power{i}(:,k) = modes(i).output_power(:,2);
        progress.out_degrees{i}(:,k) = modes(i).output_degrees;
        progress.res_norm{i}(k) = modes(i).phys_res_norm;
    end

    % Visualize epsilon.
    figure(1); title('epsilon');
    subplot 111; imagesc(modes(mode_sel(1)).epsilon{3}'); axis equal tight;

    % Visualize certain fields.
    figure(2);
    for i = mode_sel
        subplot(N, 1, i); 
        imagesc(abs(modes(i).E{3})'); 
        axis equal tight;
    end

    
    % Field design objective (power).
    figure(3); 
    my_lineplot(@plot, progress.iters, progress.out_power, mode_sel);
    title('Output powers');

    % Field design objective (angle).
    figure(4);
    my_lineplot(@plot, progress.iters, progress.out_degrees, mode_sel);
    title('Output angles');

    % Structure design objective.
    figure(5);
    my_lineplot(@semilogy, progress.iters, progress.res_norm, mode_sel);
    title('Physics residuals');

    % Physics residual.
    figure(6);
    plot(progress.struct_obj);
    title('Structure objective');

    drawnow

    progress_out = progress;
end

%% Private functions
function my_lineplot(plot_fun, x, y, mode_sel)
    subplot 111;
    for i = 1 : length(y) 
        plot_fun(x, y{i}', 'k-');
        hold on 
        if any(i == mode_sel)
            for j = 1 : size(y{i}, 1)
                plot_fun(x, y{i}(j,:), my_linestyle(i,j));
            end
        end
    end
    hold off
end


function [linestyle] = my_linestyle(i, j)
    linecolors =  'bgrcm';
    linemarkers = '.ox+*sd';
    color_ind = mod(i-1, length(linecolors)) + 1;
    marker_ind = mod(j-1, length(linemarkers)) + 1;
    linestyle = [linecolors(color_ind), linemarkers(marker_ind), '-'];
end
