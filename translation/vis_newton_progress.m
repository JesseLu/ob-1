function vis_newton_progress(progress, mode_sel)
    figure(7); 
    subplot 111;
    for i = 1 : length(progress)
        data{i} = [progress{i}(:).newton_dec];
    end
    custom_lineplot(@semilogy, [], data, mode_sel);
%     hold off
%     for i = 1 : length(progress)
%         semilogy([progress{i}(:).newton_dec], 'k.-');
%         hold on
%     end
    title('Convergence of Newton Method');
    drawnow
end


