% pop_micro_plotseg() - plots microstate segmentation over the GFP.
%
% EEGlab pop function for plotting microstate segmentation over the GFP,
% with a colour for each Microstate. Optionally the microstate numbers over
% each segment and the topographical microstate maps can by plotted.
% Calls the function MicroPlotSegments(). See this for more information. 
% 
% Usage:
%   >> pop_micro_plotseg ( EEG ); % pop up window
%   >> pop_micro_plotseg ( EEG, 'key1', 'val1', 'key2', 'val2' ... )
%
%  Please cite this toolbox as:
%  Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (unpublished
%  manuscript). Microstate EEGlab toolbox: An introductionary guide.
% 
%  Inputs
%  EEG      - EEG-lab EEG structure (channels x samples (x epochs)) with
%             the fields 'data', 'times', 'srate', and 'chanlocs'; and the 
%             'data', 'labels', and 'prototypes' fields in EEG.microstates.
%
%  Optional input:
%  'plotsegnos' - Plot Microstate numbers above microstate segments?
%                 'first' plots a number over the first segment
%                 for each microstate, 'all' (default) plots over all
%                 segments. 'none' means no numbers are plotted.
%  'plottopos'  - Plot topography of microstates? 1, yes (default); 0, no.
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

function [EEG, com] = pop_micro_plotseg(EEG, varargin)
%% Error check and initialisation
if nargin < 1
    help pop_micro_plotseg;
    return;
end;

% check whether necessary microstate substructures exist.
if ~isfield(EEG,'microstate')
   error('No microstate data present. Run microstate segmentation first.') 
end
com = '';


%% pop-up window in case no further input is given
if nargin < 2
    settings = input_popup();
    if strcmp(settings,'cancel')
        return
    end
    popup_called = 1;
else
    settings = check_settings(varargin);
    popup_called = 0;
end


%% Define command string
if popup_called
    com = sprintf('figure;MicroPlotSegments( %s', inputname(1));
else
    com = sprintf('MicroPlotSegments( %s', inputname(1));
end
com = settings_to_string(com,settings);
com = [com ' );'];


%% Run microstate fitting by evaluating com-string
eval(com)

end

function settings = input_popup()
% Function for creating popup window to input algorithm settings
%

%% Create Inputs for popup
% Label type?
style.label_type = 'popupmenu';
popmenu.label_type = {'segmentation' 'backfit'};
lab_str = popmenu.label_type{1}; %string for popupmenu
for lt = 2:length(popmenu.label_type); lab_str = [lab_str '|' popmenu.label_type{lt}]; end;
lab_tipstr = ['Should the plotted labels stem from the microstate '...
    'segmentation or from backfitting to EEG?'];
line.label_type = { {'Style' 'text' 'string' 'Plot labels from:' 'tooltipstring' lab_tipstr}, ...
    {'Style' style.label_type 'string' lab_str 'tag' 'label_type' 'value' 1} };
geo.label_type = {[1 1]};

% Plot segment numbers?
style.plotsegnos = 'popupmenu';
popmenu.plotsegnos = {'first' 'all' 'none'};
segno_str = popmenu.plotsegnos{1}; %string for popupmenu
for seg = 2:length(popmenu.plotsegnos); segno_str = [segno_str '|' popmenu.plotsegnos{seg}]; end;
line.plotsegnos = { {'Style' 'text' 'string' 'Plot microstate numbers above segments:'}, ...
    {'Style' style.plotsegnos 'string' segno_str 'tag' 'plotsegnos' 'value' 2} };
geo.plotsegnos = {[1 1]};

% Plot topographies?
style.plottopos = 'checkbox';
line.plottopos = { {'Style' style.plottopos 'value' 1 'string' 'Plot microstate topographies.' ...
    'tag' 'plottopos'} {} };
geo.plottopos = {[1 1]};


%% Order inputs for GUI
geometry = [geo.label_type geo.plotsegnos geo.plottopos];
uilist = [line.label_type line.plotsegnos line.plottopos];


%% Create Popup
[~,~,~,pop_out] = inputgui( geometry, uilist, ...
    'pophelp(''pop_micro_plotseg'');', ...
    'Draw microstate segmentation over the GFP -- pop_micro_plotseg()');


%% Interpret output from popup
if isstruct(pop_out)
    settings = struct;
    settings = interpret_popup(pop_out, settings, style, popmenu);
    
else
    settings = 'cancel';
end

end

% -------------- helper functions -------------- %
function settings = interpret_popup(pop_out, settings, style, popmenu)
% Interpret output from pop_up window, "pop_out", and arrange it in
% "settings" struct. The fields in "style" should be the same as in "pop_out"
% (defined as tags in inputgui.m). The struct popmenu is optional and only
% needed if popmenus are used in the pop_up window.

names = fieldnames(style);
for i = 1:length(names)
    switch style.(names{i})
        case 'edit'
            if isempty(pop_out.(names{i})) % empty?
                settings.(names{i}) = [];
            else
                settings.(names{i}) = eval(pop_out.(names{i}));
            end
        case 'checkbox'
            settings.(names{i}) = pop_out.(names{i});
        case 'popupmenu'
            settings.(names{i}) = popmenu.(names{i}){pop_out.(names{i})};
    end
end

end

function settings = check_settings(vargs)
%% check settings
% Checks settings given as optional inputs for MicroPlot.
% Undefined inputs is set to default values.
varg_check = {   'label_type'  'string' []  'segmentation' ;
    'plotsegnos'  'string' []  'all' ;
    'plottopos' 'integer' [] 1};
settings = finputcheck( vargs, varg_check);
if ischar(settings), error(settings); end; % check for error
end

function com = settings_to_string(com,settings)
%%
% Adds settings struct to existing com string in the form 'key1', 'val1',
% 'key2', 'val2' ... .
% Can handle structs, strings, vectors and scalars. I.e. not matrices.

names = fieldnames(settings);

for i = 1:length(names)
    if isstruct(settings.(names{i})) % struct?
        com = settings_to_string(com,settings.(names{i}));
    elseif isempty(settings.(names{i})) % empty?
        com = [ com sprintf(', ''%s'', []', names{i}) ];
    elseif ischar(settings.(names{i})) % string?
        com = [ com sprintf(', ''%s'', ''%s''', names{i}, settings.(names{i})) ];
    elseif length(settings.(names{i})) > 1 % vector?
        N_elements = length(settings.(names{i}));
        range = max(settings.(names{i})) - min(settings.(names{i})) + 1;
        if  N_elements == range % write vetor as 'min_value:max_value'
            com = [ com sprintf(', ''%s'', %g:%g', names{i}, ...
                min(settings.(names{i})), max(settings.(names{i}))) ];
        else % write vector with individual elements
            com = [ com sprintf(', ''%s'', [%g', names{i}, ...
                settings.(names{i})(1))];
            for n = 2:N_elements
                com = [ com sprintf(',%g', settings.(names{i})(n))];
            end
            com = [ com ']'];
        end
    else % scalar
        com = [ com sprintf(', ''%s'', %g', names{i}, settings.(names{i})) ];
    end
end

end