function MicroPlotFitmeas(Res, Measures, Nrange, plot_idx)
%% Plots measures of fit defined in Measures (as a cell) and contained in Res.
% If more than one measure is given, they are normalised to [0 1] and
% plotted on same plot.
%
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% July 2017.

Nmeasures = length(Measures);
hold all
for m = 1:Nmeasures
    measure = Res.(Measures{m})(plot_idx);
    
    if Nmeasures > 1
        % Normalise to [0 1]
        measure = measure - min(measure);
        measure = measure / max(measure);
    end
    
    % Plot
    plot(Nrange, measure, 'linewidth', 2)
end

xlabel('Number of Microstates')

if Nmeasures > 1
    ylabel('Normalised measure of fit (arbitrary units)')
    legend(Measures)
else
    ylabel(Measures)
end



end

