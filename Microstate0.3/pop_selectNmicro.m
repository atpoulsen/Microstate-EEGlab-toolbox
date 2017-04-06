% [] = pop_selectNmicro(EEG,epoch)
%  Select number of microstates to use. By default this function will plot
%  a new figure with four measures of fit of the best segmentation for
%  each number of microstates.
%  The user can also directly define the active number of microstates
%  without plotting measures.
%  Note: Use the function MicroPlotFitmeas to plot measures of fit in an
%  existing figure.
%
% Usage:
%   >> OUTEEG = pop_selectNmicro ( INEEG ); % pop up window
%   >> OUTEEG = pop_selectNmicro( INEEG, 'key1', 'val1', 'key2', 'val2', ....)
%
%  Note - Early untested version.
%
%  Please cite this toolbox as:
%  Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (unpublished
%  manuscript). Microstate EEGlab toolbox: An introductionary guide.
%
%  Inputs
%  EEG      - EEG-lab EEG structure (channels x samples (x epochs)) with
%             EEG ms.labels
%
%  Optional input:
%  'Nrange'      - The range of numbers of microstates to be displayed.
%  'Measures'    - Cell array of strings defining, which measures of fit to
%                  plot. Can also be a single string. Default is 'ALL',
%                  which plots all measures. Possible strings: 'KL', 'W',
%                  'CV', 'GEV' and 'ALL'.
%  'do_subplots' - If set, the created figure will contain a subplot for
%                  each measure of fit, instead of plotting all in one
%                  plot. 1 for subplots, 0 for all in on plot (default).
%  'Nmicro'      - If defined, this is made the active number of microstates
%                  used. If this is set, no figure will be created. Is
%                  empty by default. Can be an interger or empty.
%
% Outputs:
%   OUTEEG - Output dataset. This function saves output in the substruct
%            OUTEEG.microstate, which contains info general and algorithm-
%            specific settings as well as the results of segemntation.
%
% Authors:
%
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% April 2017.
%
% Copyright (C) 2017  Andreas Trier Poulsen, atpo@dtu.dk
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [OUTEEG, com] = pop_selectNmicro(EEG,varargin)
%% Error check and initialisation
com = '';

if nargin < 1
    help selectNmicro;
    return;
elseif nargin < 2
    % pop-up
    
end;
OUTEEG = EEG;
settings = check_settings(varargin, EEG);
Nmicro = settings.Nmicro;


%% Check and ordering of selected measures of fit
if ischar(settings.Measures)
    if strcmp(settings.Measures,'ALL')
        Measures = {'CV', 'GEV', 'W', 'KL'};
    else
        Measures = {settings.Measures};
    end
else
    Measures = settings.Measures;
end


% check if measures of fit are present
for m = 1:length(Measures)
    if ~isfield(EEG.microstate.Res, Measures{m})
        error('Could not find the measure %s in EEG.microstate.Res.',Measures{m})
    end
end


%% Pop-up for plotting measures of fit and selecting active number of microstates
if isempty(Nmicro)
    Nmicro = fitmeas_popup(EEG, Measures, settings);
    if strcmp(Nmicro,'cancel')
        return
    end
end


%% Selecting active number of microstates
Nmicro_ind = find(EEG.microstate.algorithm_settings.Nmicrostates == Nmicro);

if isempty(Nmicro_ind)
    error('The selected number of microstates has not been run.')
else
    OUTEEG.microstate.scalp_maps = EEG.microstate.Res.A_all{Nmicro_ind};
    OUTEEG.microstate.labels = EEG.microstate.Res.L_all{Nmicro_ind};
    OUTEEG.microstate.Res.K_act = Nmicro;
end

end

% -------------- helper functions -------------- %
function settings = check_settings(vargs, EEG)
%% check settings
% Checks settings given as optional inputs for MicroPlot.
% Undefined inputs is set to default values.
varg_check = {   'Measures'  {'string' 'cell'}  []  'ALL' ;
    'Nrange' 'real' [] EEG.microstate.algorithm_settings.Nmicrostates;
    'do_subplots' 'integer' [0 1] 0;
    'Nmicro' 'integer' [] []};
settings = finputcheck( vargs, varg_check);
if ischar(settings), error(settings); end; % check for error
end

function Nmicro = fitmeas_popup(EEG, Measures, settings)
% Creates a pop-up with plot(s) for measures of fit, and a button to
% select number of microstates.
Nmeas = length(Measures);
Nmicro = 'cancel'; % if user presses cancel button.


%% Create figure
h = figure('Units', 'normalized','position',[.2 .2 .6 .6], 'Visible','off');
if settings.do_subplots
    for m = 1:Nmeas
        subplot(Nmeas+1,1,m)
        MicroPlotFitmeas(EEG.microstate.Res, {Measures{m}}, settings.Nrange)
    end
else
    subplot('position',[0.1 .3 .85 .65])
    MicroPlotFitmeas(EEG.microstate.Res, Measures, settings.Nrange)
end


%% Nmicro input
Nmicro_txt = uicontrol('Style','text','Units', 'normalized',...
    'Position',[.2 .15 .3 .04], 'HorizontalAlignment', 'right',...
    'String','No. of microstates: ');
Nmicro_edit = uicontrol('Style','edit', 'Units', 'normalized',...
    'Position',[.51 .15 .05 .04], 'HorizontalAlignment', 'left',...
    'BackgroundColor', [1 1 1], 'String','  4  ');


%% Buttons
callback_cancel = 'close(gcf), return';
btn_cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel',...
    'Units', 'normalized', 'Position', [.69 .1 .1 .04],...
    'Callback', callback_cancel);

% callback_ok = 'handles=guidata(gcf); eval(''Nmicro = get( handles.Nmicro_edit, ''string'')''); close(gcf)';
% callback_ok = 'close(gcf)';
callback_ok = 'set(gcbo, ''userdata'', 1);';
btn_ok = uicontrol('Style', 'pushbutton', 'String', 'OK',...
    'Units', 'normalized', 'Position', [.8 .1 .1 .04],...
    'tag', 'ok', 'Callback', callback_ok);

set(h,'Visible','on');
% guidata(h,handles); % to be able to access Nmicro_edit, when callin from within a function
waitfor( findobj('parent', h, 'tag', 'ok'), 'userdata');
Nmicro = eval(get(Nmicro_edit, 'string'));


close(h);


end
