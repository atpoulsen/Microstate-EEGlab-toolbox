% pop_micro_import_proto() - import maps from other dataset.
%
% Usage:
%   >> EEG = pop_micro_import_proto ( EEG, ALLEEG ); % pop up window
%   >> EEG = pop_micro_import_proto ( EEG, ALLEEG, data_idx )
%
% Please cite this toolbox as:
% Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (unpublished
% manuscript). Microstate EEGlab toolbox: An introductionary guide.
%
% Inputs:
%   EEG             - EEGlab EEG structure.
%   ALLEEG          - EEGlab structure containing all datasets read into
%                     EEGlab as EEG structures.
%
% Optional inputs:
%  dataset_idx      - Index to select which dataset in AALEEG to import
%                     microstate maps from.
%
% Outputs:
%   EEG             - EEG-lab EEG structure containing new microstate
%                     prototypes in EEG.microstate.prototypes.
%
% Authors:
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% May 2017.
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
function [EEG, com] = pop_micro_import_proto(EEG, ALLEEG, dataset_idx)
%% Error check and initialisation
com = '';

if nargin < 2
    help pop_micro_import_proto
    return;
end


%% Pop-up
if nargin < 3
    % pop-up window in case no further input is given
    dataset_names = {ALLEEG.setname};
    [dataset_idx, confirmed] = listdlg2('ListString',dataset_names,'PromptString',...
        {'Select dataset to import prototypes from:'}, 'selectionmode', 'single');
    if ~confirmed
        return
    end
end


%% Import microstate maps
if isfield(EEG,'microstate')
    if isfield(EEG.microstate,'prototypes')
        warning('Overwriting existing prototypes in dataset.')
    end
end

if ~isfield(ALLEEG(dataset_idx).microstate,'prototypes')
   error('prototypes not present in dataset chosen to import from.') 
end

EEG.microstate.prototypes = ALLEEG(dataset_idx).microstate.prototypes;

if ~isfield(EEG.microstate,'data')
    % To avoid having to run pop_micro_selectdata 
    EEG.microstate.data = 'set_data';
end

%% Define command string
com = sprintf('%s = pop_micro_import_proto( %s, %s, %d)', inputname(1), ...
    inputname(1), inputname(2), dataset_idx);

end