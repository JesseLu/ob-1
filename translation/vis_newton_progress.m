function vis_newton_progress(progress, mode_sel)
    subplot(2, 6, 6+5);
    for i = 1 : length(progress)
        data{i} = [progress{i}(:).newton_dec];
    end
    custom_lineplot(@semilogy, [], data, mode_sel);
    title('Convergence of Newton Method');
    drawnow

    subplot(2, 6, 6+6);
end


