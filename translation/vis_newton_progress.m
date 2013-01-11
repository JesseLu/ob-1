function vis_newton_progress(progress)
    figure(7); 
    subplot 111;
    hold off
    for i = 1 : length(progress)
        semilogy([progress{i}(:).newton_dec], 'k.-');
        hold on
    end
    title('Convergence of Newton Method');
    drawnow
end


