% pop_micro_fit() - backfits microstate prototype maps to EEG data.
%  
% Backfit microstate prototype maps to EEG stored in EEG.data. 
%
% Usage:
%   >> EEG = pop_micro_fit ( EEG ); % pop up window
%   >> EEG = pop_micro_fit ( EEG, 'key1', 'val1')
%
% Please cite this toolbox as:
% Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (unpublished
% manuscript). Microstate EEGlab toolbox: An introductionary guide.
%
% Inputs:
%   EEG - EEGlab EEG structure with microstate field.
%
% Optional input:
%  'polarity' - Account for polarity when fitting Typically off for
%               spontaneous EEG and on for ERP data (default = 0).
%
% Outputs:
%   EEG.microstate.fit - Struct in EEG structure containing info on
%                        microstate fitting:
%                     .labels   - Microstate labels (N microstates x time
%                                 (x trials)).
%                     .polarity - See inputs.
%
% Authors:
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% Andreas Pedroni, andreas.pedroni@uzh.ch
% University of Zürich, Psychologisches Institut, Methoden der
% Plastizitätsforschung. 
%
% September 2017.
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

function [EEG, com] = pop_micro_fit(EEG, varargin)
%% Error check and initialisation
if nargin < 1
    help pop_micro_fit;
    return;
end

% check whether microstate substruct exists.
if ~isfield(EEG,'microstate')
   error('No microstate data present. Run microstate segmentation first.') 
end

com = '';


%% pop-up window in case no further input is given
if nargin < 2
    settings = fit_popup();
    if strcmp(settings,'cancel')
        return
    end
    disp('Fitting microstate prototype maps to EEG data...')
else
    settings = check_settings(varargin);
end


%% Run backfitting
L = MicroFit(EEG.data, EEG.microstate.prototypes, settings.polarity);
EEG.microstate.fit.labels = L;
EEG.microstate.fit.polarity = settings.polarity;


%% Define command string
com = sprintf('%s = pop_micro_fit( %s', inputname(1), inputname(1));
com = settings_to_string(com,settings);
com = [com ' );'];

end

% ------------------------------ Pop-up ---------------------------------%
function settings = fit_popup()
% Function for creating popup window to input settings
%

%% Create Inputs for popup
% Polarity
% style.polarity = 'checkbox';
% pol_tipstr = 'Spontaneous EEG typically ignore polarity (off). Typically on for ERP data';
% line.polarity = { {'Style' style.polarity 'value' 0 'string' 'Account for polarity when fitting ' ...
%     'tooltipstring' pol_tipstr 'tag' 'polarity'} {} };
% geo.polarity = {[1 1]};
style.polarity = 'checkbox';
pol_tipstr = 'Spontaneous EEG typically ignore polarity (off). Typically on for ERP data.';
line.polarity = { {'Style' 'text' 'string' 'Account for polarity when fitting.' ...
    'tooltipstring' pol_tipstr}, ...
    {'Style' style.polarity 'value' 0 'tag' 'polarity'} };
geo.polarity = {[1 .2]};


%% Order inputs for GUI
geometry = geo.polarity;
uilist = line.polarity;


%% Create Popup
[~,~,~,pop_out] = inputgui( geometry, uilist, ...
    'pophelp(''pop_micro_fit'');', 'Fit microstate maps to EEG -- pop_micro_fit()');
 

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
function settings = check_settings(vargs)
% Checks settings given as optional inputs for MicroFit.
% Undefined inputs is set to default values.
%%
varg_check = { 'polarity'      'integer'	[]	0 };

settings = finputcheck( vargs, varg_check);
if ischar(settings), error(settings); end % check for error
end

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

function com = settings_to_string(com,settings)
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