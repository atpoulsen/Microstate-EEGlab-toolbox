% pop_micro_smooth() - Temporally smooths microstate labels
%  
% Temporally smooths microstate labels obtained either from segmentation or
% from backfitting. The labels should be stored in either 
% EEG.microstate.labels or EEG.microstate.fit.labels, correspondingly.
%
% Usage:
%   >> EEG = pop_micro_smooth ( EEG ); % pop up window
%   >> EEG = pop_micro_smooth ( EEG, 'key1', 'val1')
%
% Please cite this toolbox as:
% Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (unpublished
% manuscript). Microstate EEGlab toolbox: An introductionary guide.
%
% Inputs:
%   EEG  - EEGlab EEG structure with microstate field.
%
% Optional input:
%  'label_type'  - Smooth labels obtained from: 'segmentation' (default) 
%                  or 'backfit'.
%  'smooth_type' - Smoothing type: 'reject segments' (default) or
%                  'windowed'.
%
% Method-specific inputs:
% * Reject segments:
%   'minTime'  - Redristibute segments smaller than minTime (in ms) to
%                the next best fitting microstate (default = 20 ms).
%   'polarity' - Account for polarity when calculating the global map
%                dissimilarity. Typically off for spontaneous EEG and on
%                for ERP data (default = 0).
% * Windowed smoothing:
%   'smooth_width'   - Temporal smoothing width. Integer denoting the 
%                      number of samples on each side of current sample
%                      (default: 3).
%   'smooth_weight'  - Temporal smoothing weight (default: 5).
%   'max_iterations' - Maximum number of iterations of algorithm
%                      (default: 1000).
%   'threshold'      - Threshold of convergence based on relative change
%                      in noise variance (default: 1e-6).
%
% Output:
%   EEG - EEGlab struct with updated microstate labels.
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
% See also: eeglab pop_micro_segment pop_micro_fit

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

function [EEG, com] = pop_micro_smooth(EEG, varargin)
%% Error check and initialisation
if nargin < 1
    help pop_micro_smooth;
    return;
end

% check whether microstate substruct exists.
if ~isfield(EEG,'microstate')
   error('No microstate data present. Run microstate segmentation first.') 
end

com = '';


%% pop-up window in case no further input is given
if nargin < 2
    settings = smooth_popup();
    if strcmp(settings,'cancel')
        return
    end
    if strcmp(settings.smooth_type,'reject segments')
        disp('Smoothing microstate labels using reject small segments...')
    else
        disp('Smoothing microstate labels using windowed smoothing...')
    end
else
    settings = check_settings(varargin);
end


%% Select data
if strcmp(settings.label_type,'segmentation')
    if ischar(EEG.microstate.data)
        data = EEG.data;
    else
        data = EEG.microstate.data;
    end
else
    data = EEG.data;
end
    
% Prototypes
A = EEG.microstate.prototypes;

%% Ready smoothing-specific settings
if strcmp(settings.smooth_type,'windowed')
    opts.b = settings.smooth_width;
    opts.lambda = settings.smooth_weight;
    opts.max_iterations = settings.max_iterations;
    opts.thresh = settings.threshold;
else
    opts.polarity = settings.polarity;
    
    % converting minTime from ms to samples
    minTime_ms = settings.minTime;
    minTime_samples = round( minTime_ms * EEG.srate/1000 );
    opts.minTime = minTime_samples;
end


%% Run smoothing
labels = MicroSmooth(data, A, settings.smooth_type, opts);


%% Save smoothed labels
if strcmp(settings.label_type,'segmentation')
    EEG.microstate.labels = labels;
else
    EEG.microstate.fit.labels = labels;
end


%% Define command string
com = sprintf('%s = pop_micro_smooth( %s', inputname(1), inputname(1));
com = settings_to_string(com,settings);
com = [com ' );'];
end

% ------------------------------ Pop-ups ---------------------------------%
function settings = smooth_popup()
% Function for creating popup window to input settings
%

%% Create Inputs for popup
% Select labels to smooth
style.label_type = 'popupmenu';
dropdown_label = {'Microstate segmentation' 'Backfitting prototypes to EEG'}; % For dropdown menu
popmenu.label_type = {'segmentation' 'backfit'}; % Corresponding calls for pop-function
label_str = dropdown_label{1}; %string for popupmenu
for l = 2:length(dropdown_label); label_str = [label_str '|' dropdown_label{l}]; end
line.label_type = { {'Style' 'text' 'string' 'Smooth labels obtained from:'}, ...
    {'Style' style.label_type 'string' label_str 'tag' 'label_type' 'value' 1} };
geo.label_type = {[1 1]};

% Smoothing type
style.smooth_type = 'popupmenu';
dropdown_smooth = {'Reject small segments' 'Windowed smoothing'}; % For dropdown menu
popmenu.smooth_type = {'reject segments' 'windowed'}; % Corresponding calls for pop-function
smooth_str = dropdown_smooth{1}; %string for popupmenu
for s = 2:length(dropdown_smooth); smooth_str = [smooth_str '|' dropdown_smooth{s}]; end
line.smooth_type = { {'Style' 'text' 'string' 'Choose smoothing method:'}, ...
    {'Style' style.smooth_type 'string' smooth_str 'tag' 'smooth_type' 'value' 1} };
geo.smooth_type = {[1 1]};


%% Order inputs for GUI
geometry = [geo.label_type geo.smooth_type];
uilist = [line.label_type line.smooth_type];


%% Create Popup
[~,~,~,pop_out] = inputgui( geometry, uilist, ...
    'pophelp(''pop_micro_smooth'');', 'Smooths microstate labels -- pop_micro_smooth()');
 

%% Interpret output from popup
if isstruct(pop_out)
    settings = struct;
    settings = interpret_popup(pop_out, settings, style, popmenu);
else
    settings = 'cancel';
end


%% Create popup for smoothing specific settings
if ~strcmp(settings,'cancel')
    switch settings.smooth_type
        case 'reject segments'
            settings = reject_popup(settings);
        case 'windowed'
            settings = window_popup(settings);
        otherwise
            error(['selected smooth type,''' settings.smooth_type ''' not available'])
    end
end
end

function settings = reject_popup(settings)
% Popup for unique input for reject small segments
%

%% Create Inputs for popup
% Title string
info_str = 'Input parameters specific for ''Reject small segments''.';
line.info = { {'Style' 'text' 'string' info_str} {} };
geo.info = {1 1};

% Redristibute segments (minTime)
style.minTime = 'edit';
line.minTime = { {'Style' 'text' 'string' 'Redistribute segments smaller than (in ms):'}, ...
    {'Style' style.minTime 'string' ' 30 ' 'tag' 'minTime'} };
geo.minTime = {[1 .2]};

% Polarity
style.polarity = 'checkbox';
pol_tipstr = 'Spontaneous EEG typically ignore polarity (off). Typically on for ERP data.';
line.polarity = { {'Style' 'text' 'string' 'Account for polarity when fitting.' ...
    'tooltipstring' pol_tipstr}, ...
    {'Style' style.polarity 'value' 0 'tag' 'polarity'} };
geo.polarity = {[1 .2]};


%% Order inputs for GUI
geometry = [geo.info geo.minTime geo.polarity];
uilist = [line.info line.minTime line.polarity];


%% Create Popup
[~,~,~,pop_out] = inputgui( geometry, uilist, ...
    'pophelp(''pop_micro_smooth'');', 'Extra input for reject small segments -- pop_micro_smooth()');


%% Interpret output from popup
if isstruct(pop_out)
    settings = interpret_popup(pop_out, settings, style);
else
    settings = 'cancel';
end
end

function settings = window_popup(settings)
% Popup for unique input for windowed smoothing
%

%% Create Inputs for popup
% Title string
info_str = 'Input parameters specific for ''Windowed smoothing''.';
line.info = { {'Style' 'text' 'string' info_str} {} };
geo.info = {1 1};

% Smoothing width
style.smooth_width = 'edit';
line.smooth_width = { {'Style' 'text' 'string' 'Width of smoothing windows (in samples):'}, ...
    {'Style' style.smooth_width 'string' ' 3 ' 'tag' 'smooth_width'} };
geo.smooth_width = {[1 .2]};

% Smoothing weight
style.smooth_weight = 'edit';
line.smooth_weight = { {'Style' 'text' 'string' 'Smoothing weight:'}, ...
    {'Style' style.smooth_weight 'string' ' 5 ' 'tag' 'smooth_weight'} };
geo.smooth_weight = {[1 .2]};

% Max iterations
style.max_iterations = 'edit';
line.max_iterations = { {'Style' 'text' 'string' 'Max. no. of iterations:'}, ...
    {'Style' style.max_iterations 'string' ' 1000 ' 'tag' 'max_iterations'} };
geo.max_iterations = {[1 .2]};

% Threshold
style.threshold = 'edit';
line.threshold = { {'Style' 'text' 'string' 'Relative threshold for convergence:'}, ...
    {'Style' style.threshold 'string' ' 1e-6 ' 'tag' 'threshold'} };
geo.threshold = {[1 .2]};


%% Order inputs for GUI
geometry = [geo.info geo.smooth_width geo.smooth_weight geo.max_iterations geo.threshold];
uilist = [line.info line.smooth_width line.smooth_weight line.max_iterations line.threshold];


%% Create Popup
[~,~,~,pop_out] = inputgui( geometry, uilist, ...
    'pophelp(''pop_micro_smooth'');', 'Extra input for windowed smoothing -- pop_micro_smooth()');


%% Interpret output from popup
if isstruct(pop_out)
    settings = interpret_popup(pop_out, settings, style);
else
    settings = 'cancel';
end
end
% ----------------------------------------------------------------------- %

% -------------------------- Helper functions --------------------------- %
function settings = check_settings(vargs)
% Checks settings given as optional inputs for pop_micro_smooth.
% The function checks and rearranges optional inputs to pop_micro_smooth.m
% struct. Undefined inputs is set to default values.
%% General inputs
varg_check = { 'label_type'  'string'    []         'segmentation'
    'smooth_type'  'string'    []         'reject segments'};


%% Algorithm specifik inputs
smooth_type = find(strcmp('smooth_type',vargs)) + 1;
switch vargs{smooth_type}
    case 'reject segments'
        varg_check = [varg_check;
            {'minTime'  'integer'    []         30;
            'polarity'  'integer'    []         0 } ];
    case 'windowed'
        varg_check = [varg_check;
            { 'smooth_weight'  'integer'    []         5;
            'smooth_width'  'integer'    []         3;
            'max_iterations'  'integer'    []         1000;
            'threshold'  'real'    []         1e-6 } ];
    otherwise
        error(['selected algorithm,''' vargs{smooth_type} ''' not available'])
end
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