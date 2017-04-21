% eegplugin_microstate()
%
% Usage:
%   >> eegplugin_microstate(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Please cite this toolbox as:
% Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (unpublished
% manuscript). Microstate EEGlab toolbox: An introductionary guide.
%
% Authors:
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% Andreas Pedroni, andreas.pedroni@uzh.ch
% University of Zürich, Psychologisches Institut, Methoden der
% Plastizitätsforschung. 
%
% April 2017.
%
% See also: eeglab

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

function eegplugin_microstate(fig, trystrs, catchstrs)
if nargin < 3
    error('eegplugin_microstate requires 3 arguments');
end

%% Create menus
% Tools menu
toolsmenu = findobj(fig, 'tag', 'tools');
submenutools = uimenu( toolsmenu, 'label', 'Microstate analysis', 'separator', 'on');

% % Plot menu
% plotmenu = findobj(fig, 'tag', 'plot');
% submenuplot = uimenu( plotmenu, 'label', 'Microstate analysis', 'separator', 'on');


%% Menu callback commands for submenu
% Analysis functions
complaceholder = 'error(''sorry function not implemented yet. This is just a placeholder'');';
comdata = [trystrs.no_check complaceholder catchstrs.store_and_hist];
comsegment_simple = [trystrs.no_check complaceholder catchstrs.store_and_hist];
comsegment_advanced = [trystrs.no_check '[EEG LASTCOM]=pop_micro_segment(EEG);' ...
    catchstrs.store_and_hist];
combackfit = [trystrs.no_check '[EEG LASTCOM]=pop_micro_fit(EEG);' ...
    catchstrs.store_and_hist];
comselect = [trystrs.no_check '[EEG LASTCOM]=pop_selectNmicro(EEG);' ...
    catchstrs.store_and_hist];
comstats = [trystrs.no_check '[EEG LASTCOM]=pop_micro_stats(EEG);' ...
    catchstrs.store_and_hist];

% Plot functions
scalp_str = 'LASTCOM = [''figure;topo_micro(EEG.microstate.scalp_maps,EEG.chanlocs);'']; eval(LASTCOM);';
complotscalp = [trystrs.no_check scalp_str catchstrs.store_and_hist];
complotseg = [trystrs.no_check '[EEG,LASTCOM] = pop_micro_plotseg(EEG);' ...
    catchstrs.store_and_hist];


%% Create Tools submenus
% Simple segmentation options
uimenu( submenutools, 'Label', 'Segment into microstates (Quick start) (placeholder)',...
    'CallBack', comsegment_simple);

% Advanced segmentation options
uimenu( submenutools, 'Label', 'Select data for analysis (advanced settings) (placeholder)', ...
    'CallBack', comdata,'separator', 'on');
uimenu( submenutools, 'Label', 'Segment into microstates (advanced settings)', ...
    'CallBack', comsegment_advanced);

% Post segmentation
uimenu( submenutools, 'Label', 'Select active number of microstates', ...
    'CallBack', comselect,'separator', 'on');
uimenu( submenutools, 'Label', 'Backfit microstates on EEG', ...
    'CallBack', combackfit);
uimenu( submenutools, 'Label', 'Calculate microstate statistics', ...
    'CallBack', comstats);

% Plot functions
uimenu( submenutools, 'Label', 'Plot microstates scalp topographies', 'CallBack',...
    complotscalp,'separator', 'on');
uimenu( submenutools, 'Label', 'Plot microstate segmentations','CallBack', complotseg);


% %% Create Plot submenus
% uimenu( submenuplot, 'Label', 'Plot microstates scalp topographies', 'CallBack', complotscalp);
% uimenu( submenuplot, 'Label', 'Plot microstate segmentation','CallBack', complotseg);

end