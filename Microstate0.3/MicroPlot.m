%   [] = MicroPlot(EEG,epoch)
%   Draws a figure with Microstate segments over the GFP
%   >> OUTEEG = MicroPlot( INEEG, 'key1', 'val1', 'key2', 'val2', ....)
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
%  'epoch'     - timewindow of analysis (vector of timeframes).
%  'Nplots'    - number of subplots. Has to be <= size(EEG.data,3)
%
% Authors:
%
% Andreas Pedroni, andreas.pedroni@uzh.ch
% University of Zürich, Psychologisches Institut, Methoden der
% Plastizitätsforschung.
%
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% February 2017.

% Copyright (C) 2017 Andreas Pedroni, andreas.pedroni@uzh.ch.
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

function MicroPlot(EEG,varargin)

%% Error check and initialisation
if nargin < 1
    help MicroPlot;
    return;
end;

settings = check_settings(varargin, EEG);


%%
GFP = EEG.microstate.fit.GFP(:,settings.epoch);
Label = EEG.microstate.fit.bestLabel(:,settings.epoch);

   
c20 =   [0.368627450980392,0.309803921568627,0.635294117647059;0.232941976907398,0.426964379443341,0.716724269529944;
        0.197103273996820,0.545428702089154,0.740441532379098;0.289164658427885,0.676519464158588,0.684127372366063;
        0.425488941852069,0.775245691220066,0.646374895135081;0.570745493485241,0.831883120082506,0.644760405293702;
        0.712331911655248,0.884123184534485,0.639177467615844;0.852572369568295,0.945354582576408,0.606699084357464;
        0.916245263768755,0.960207509159080,0.609279235239802;0.939109517434401,0.954525564903514,0.691784259410135;
        0.964614344304082,0.940412515225342,0.691500293230307;0.992457455354208,0.899715439687825,0.585852977458640;
        0.995633845898762,0.831741829511967,0.491082885642943;0.993209665846219,0.718091384842371,0.403306245777860;
        0.984967339319457,0.589815920663289,0.324081427293180;0.962493646532822,0.451160196578616,0.265020625541367;
        0.918594329250906,0.348117017024907,0.280748065372818;0.843470509188835,0.253888856898642,0.309426572786864;
        0.746582251164705,0.137135217599700,0.298133552881715;0.619607843137255,0.00392156862745098,0.258823529411765];

if  size(EEG.microstate.fit.MStemplate,2) < 9
    colors = colormap('lines');
else
    colors = c20;
end

%% Plotting
% creating time axis
Nsamples = size(GFP,2);
time_axis = 1:Nsamples;
if isfield(EEG,'srate')
   time_axis = time_axis / EEG.srate; 
end

for t=1:settings.Nplots
    subplot(settings.Nplots,1,t)
    hold on
    for i=unique(Label(t,:))
        x = nan(1,Nsamples);
        idx = Label(t,:) == i; % finding indices where microstate is active
        % adding one sample to start of each segment to avoid white holes
        % between segments.
        idx = [idx(2:end) false] | idx; 

        x(idx) = GFP(t,idx);
        a = area(time_axis,x);
        a.EdgeColor = colors(i,:);
        a.FaceColor = colors(i,:);
    end
    xlim([min(time_axis) max(time_axis)])
    xlabel('Time (s)')
    ylabel('GFP')
end

end

% -------------- helper functions -------------- %
function settings = check_settings(vargs, EEG)
%% check settings
% Checks settings given as optional inputs for MicroPlot.
% Undefined inputs is set to default values.
varg_check = {   'epoch'  'integer'    []         1:size(EEG.data,2) ;
    'Nplots' 'integer' [] size(EEG.microstate.fit.bestLabel,1)};
settings = finputcheck( vargs, varg_check);
if ischar(settings), error(settings); end; % check for error
end
