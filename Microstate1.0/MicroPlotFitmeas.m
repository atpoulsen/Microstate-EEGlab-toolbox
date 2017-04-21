function MicroPlotFitmeas(Res, Measures, Nrange)
%% Plots measures of fit defined in Measures (as a cell) and contained in Res.
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% April 2017.

hold all
for m = 1:length(Measures)   
    measure = Res.(Measures{m});
    % Normalise to [0 1]
    measure = measure - min(measure);
    measure = measure / max(measure);
    
    % Plot
    plot(Nrange, measure, 'linewidth', 2)
end

xlabel('Number of Microstates')
ylabel('Normalised measure of fit (arbitrary units)')
legend(Measures)



end

