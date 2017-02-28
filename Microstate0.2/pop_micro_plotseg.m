% pop_micro_plotseg() - plots Microstate segments over the GFP.
%
%  Note - Early untested version.
%
% Usage:
%   >> pop_micro_plotseg ( EEG ); % pop up window
%   >> pop_micro_plotseg ( EEG, 'key1', 'val1', 'key2', 'val2' ... )
%
% Please cite this toolbox as:
% Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (unpublished
% manuscript). Microstate EEGlab toolbox: An introductionary guide.
%
% Inputs:
%   EEG        - EEG-lab EEG structure (channels x samples (x epochs)) with
%                .microstate.fit.bestLabel (created by MicroFit.m)
%
% Optional inputs:
%  'epoch'     - timewindow of analysis (vector of timeframes).
%  'Nplots'    - number of subplots. Has to be <= size(EEG.data,3)
%
% Authors:
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% Andreas Pedroni, andreas.pedroni@uzh.ch
% University of Zürich, Psychologisches Institut, Methoden der
% Plastizitätsforschung. 
%
% February 2017.
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
if ~isfield(EEG.microstate,'fit')
   error('No microstate fitting info present. Run microstate fit first.') 
end

com = '';


%% pop-up window in case no further input is given
if nargin < 2
    settings = plotseg_popup();
    if strcmp(settings,'cancel')
        return
    end
else
    settings = check_settings(varargin,EEG);
end

%% Remove empty fields in settings
names = fieldnames(settings);
for n = 1:length(names)
    if isempty(settings.(names{n}))
        % removing epoch field if empty
        settings = rmfield(settings,names{n});
    end
end


%% Define command string
com = sprintf('MicroPlot( %s', inputname(1));
com = settings_to_string(com,settings);
com = [com ' );'];


%% Run microstate fitting by evaluating com-string
disp('Plotting microstate segments over the GFP...')
eval(com)

end

% ------------------------------ Pop-up ---------------------------------%
function settings = plotseg_popup()
% Function for creating popup window to input settings
%

%% Create Inputs for popup
% Title string
info_str1 = 'Please note that this is an early version of the plugin. Bug-reports and suggestions';
info_str2 = 'are welcome at atpo@dtu.dk.';
line.info = { {'Style' 'text' 'string' info_str1} ...
    {'Style' 'text' 'string' info_str2} {} };
geo.info = {1 1 1};


% Epoch
style.epoch = 'edit';
epoch_tipstr = 'Vector of timeframes. Leave empty to use entire range.';
line.epoch = { {'Style' 'text' 'string' 'Timewindow of analysis', ...
    'tooltipstring' epoch_tipstr}, ...
    {'Style' style.epoch 'string' '' 'tag' 'epoch'} };
geo.epoch = {[1 .2]};

% Nplots
style.Nplots = 'edit';
plot_tipstr = 'Has to be <= size(EEG.data,3). Leave empty to set to size(EEG.data,3).';
line.Nplots = { {'Style' 'text' 'string' 'Number of subplots', ...
    'tooltipstring' plot_tipstr}, ...
    {'Style' style.Nplots 'string' '' 'tag' 'Nplots'} };
geo.Nplots = {[1 .2]};


%% Order inputs for GUI
geometry = [geo.info geo.epoch geo.Nplots];
uilist = [line.info line.epoch line.Nplots];


%% Create Popup
[~,~,~,pop_out] = inputgui( geometry, uilist, 'pophelp(''pop_micro_plotseg'');',...
    'plots microstate segments over the GFP -- pop_micro_plotseg()');
 

%% Interpret output from popup
if isstruct(pop_out)
    settings = struct;
    settings = interpret_popup(pop_out, settings, style);
else
    settings = 'cancel';
end
end
% ----------------------------------------------------------------------- %

% -------------------------- Helper functions --------------------------- %
function settings = check_settings(vargs, EEG)
%% check settings
% Checks settings given as optional inputs for MicroPlot.
% Undefined inputs is set to default values.
varg_check = {   'epoch'  'integer'    []         1:size(EEG.data,2) ;
    'Nplots' 'integer' [] size(EEG.microstate.fit.bestLabel,1)};
settings = finputcheck( vargs, varg_check);
if ischar(settings), error(settings); end; % check for error
end

function settings = interpret_popup(pop_out, settings, style, popmenu)
%%
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