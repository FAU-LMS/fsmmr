%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2021-2022 Chair of Multimedia Communications and Signal Processing (LMS),
% Friedrich-Alexander-University Erlangen-Nürnberg (FAU).
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its contributors
%    may be used to endorse or promote products derived from this software without
%    specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
% BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% This is the main.m file for Color Super Resolution of Point Clouds using
% Frequency-Selective Mesh-to-Mesh Resampling.
% Application: Upsampling in 3D space of
% point clouds with following quality evaluation
% Author: Viktoria Heimann


close all
clear

% Set parameters
interpolationMethods = {'linear', 'cubic', 'natural', 'FSMMR'}; % FSMMR from the MMSP'21 publication
cloud_name = 'Asterix';

% Initialisation of Parameters if upsampling
sample_percent = 50;
sample_ratio = sample_percent/100;

% Initialisation of FSMMR Parameters
params.blkSize = 4;
params.borderWidth = 0; % for point clouds borderWidth has to be set to 0
params.rho = 1;
params.sigma = 0.5;
params.max_iter = 1000;
params.orthogonality_correction = 0.5;
params.transform_size = 8;
params.transform_type = 'dct';

%%

% Load point cloud
pcloud_orig = pcread(['data/' cloud_name '.ply']);

% Downsample point cloud
pcloud_down = pcdownsample(pcloud_orig, 'random', sample_ratio);

% Extract point cloud points that have to be estimated for upsampled
% version
loc_orig = pcloud_orig.Location;
loc_down = pcloud_down.Location;
col_down = pcloud_down.Color;

idx_2upsample = ismember(loc_orig, loc_down, 'rows'); % all points that are not members of original and downsampled pcloud have to be upsampled
loc_2upsample = loc_orig(~idx_2upsample, :);
col_2upsample = zeros(size(loc_2upsample, 1), size(col_down, 2));
% create pcloud
pcloud_2upsample = pointCloud(loc_2upsample, 'Color', col_2upsample/255); %Contains locations for whom color should be estmated = grid

%% Here loop through interpolationMethods with same input has to be established
for ii = 1 : length(interpolationMethods)
    %% Determine interpolation method
    interpolationMethod = interpolationMethods{ii};
    
    %% Reconstruction
    pcloud_final = pcloud_upsample_color(pcloud_down, pcloud_2upsample,...
        params, interpolationMethod);    
    
    %% Save
    % define save path
    save_path = ['results/' interpolationMethod '/'];
    if ~exist(save_path, 'dir')
        mkdir(save_path)
    end
    % write
    pcwrite(pcloud_final, [save_path '/' cloud_name '.ply'])
    
    %% Show
    figure
    pcshow(pcloud_final, 'MarkerSize', 50)
    title(['2D ' interpolationMethod])
end % interpolationMethod
