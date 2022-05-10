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

%% Function to reconstruct the color attribute of point cloud objects
% in order to re- and upsample point cloud objects

function pcloud_return = reconstruct_pcloud_upsample(blk_in, ...
    processingOrder, params, blk_2upsample, interpolationMethod)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% blk_in: point cloud encapsulated in blocks that contains original
% information
% processingOrder: matrix that contains order for processing blocks
% params: containing FSMMR parameters
% blk_2upsample: point cloud encapsulated in blocks that contains locations
% for whom color values are searched.
%
% pcloud_return: point cloud on locations from blk_2upsample and newly
% generated according color values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

transform_size = params.transform_size;
transform_type = params.transform_type;
sigma = params.sigma;
blkSize = params.blkSize;
max_iter = params.max_iter;
orthogonality_correction = params.orthogonality_correction;

pcloud_return = cell(size(blk_in));
blk_return = cell(size(blk_in));

numBlocks = length(processingOrder(1,:));

for blkIterator = 1 : numBlocks % iterate through all blocks
    
    % Get next block to process, get original information of downsampled
    % point cloud
    index = [processingOrder(1,blkIterator) processingOrder(2,blkIterator)...
        processingOrder(3, blkIterator)]; %xIndex, yIndex, zIndex
    % Get information about the original points incl support area
    blk_tmp = blk_in{index(1), index(2), index(3)}; % x, y, z
    blk_2upsample_tmp = blk_2upsample{index(1), index(2), index(3)}; %x,y,z
    % Display
    if mod(blkIterator, round(numBlocks/10)) == 0
        disp(['Done ..... ' num2str(round(100*blkIterator/numBlocks)) '%'])
    end
    if isempty(blk_tmp)
        continue
    end
    if ~strcmp(interpolationMethod, 'FSMMR')
        if size(blk_tmp, 1) < 3
            continue
        end
    end
    if isempty(blk_2upsample_tmp)
        continue
    end
    % Get current blk information from downsampled pcloud
    x_tmp = blk_tmp(:,1);
    y_tmp = blk_tmp(:,2);
    z_tmp = blk_tmp(:, 3);
    RGB_tmp = blk_tmp(:, 4:end-1);
    w_tmp = blk_tmp(:, end);
    
    % Get current blk information from final pcloud
    x_2upsample_tmp = blk_2upsample_tmp(:,1);
    y_2upsample_tmp = blk_2upsample_tmp(:,2);
    z_2upsample_tmp = blk_2upsample_tmp(:,3);
    
    % Define which plane should be reconstructed
    % The plane with the smallest variance should be reconstructed, hence,
    % the other two axis are used as fixed values
    x_var = var(x_tmp);
    y_var = var(y_tmp);
    z_var = var(z_tmp);
    
    [~, idx_var] = min([x_var, y_var, z_var]);
    if idx_var == 1
        dim1 = y_tmp;
        dim2 = z_tmp;
        dim3 = x_tmp;
        
        dim1_2upsample = y_2upsample_tmp;
        dim2_2upsample = z_2upsample_tmp;
        dim3_2upsample = x_2upsample_tmp;
    elseif idx_var == 2
        dim1 = x_tmp;
        dim2 = z_tmp;
        dim3 = y_tmp;
        
        dim1_2upsample = x_2upsample_tmp;
        dim2_2upsample = z_2upsample_tmp;
        dim3_2upsample = y_2upsample_tmp;
    elseif idx_var == 3
        dim1 = x_tmp;
        dim2 = y_tmp;
        dim3 = z_tmp;
        
        dim1_2upsample = x_2upsample_tmp;
        dim2_2upsample = y_2upsample_tmp;
        dim3_2upsample = z_2upsample_tmp;
    end
    
    
    %% COLOR %%
    dim_out = [dim1_2upsample(:), dim2_2upsample(:), dim3_2upsample(:)];
    dim_in = [dim1(:), dim2(:), dim3(:)];
    
    [points2_in, points2_out] = ...
        transformSurface2Plane_shortestWay(dim_in, dim_out, transform_size, RGB_tmp);
    
    %% FSMMR for color resampling
    % Calculate basis functions for (original) points
    [basis_functions, freq_weights] = generate_floating_basis_functions...
        (points2_in(:,1), points2_in(:,2), transform_size, transform_type, 5, sigma);
    
    % Calculate basis functions for resampled points
    [basis_functions_out, ~] = generate_floating_basis_functions...
        (points2_out(:,1), points2_out(:,2), transform_size, transform_type, 1, sigma);
    
    % Apply FSMMR for color resampling
    for col_channel = 1 : 3
        if strcmp(interpolationMethod, 'FSMMR')
            used_freqs = fsmmr(RGB_tmp(:,col_channel), w_tmp,  max_iter, ...
                orthogonality_correction, basis_functions, freq_weights);
            % Use estimated expansion coefficients to evaluate on voxelized points
            points2_out(:,2+col_channel) = uint8(basis_functions_out * used_freqs(:));
        else
            points2_out(:, 2+col_channel) = griddata(double(points2_in(:, 1)), double(points2_in(:, 2)),...
                double(RGB_tmp(:, col_channel)), double(points2_out(:,1)), double(points2_out(:,2)), ...
                interpolationMethod);
        end
    end
    %     Reverse change in dimensions
    if idx_var == 1
        blk_return{index(1), index(2), index(3)} = ...
            [(index(1)-1)*blkSize + dim3_2upsample(:), ...
            (index(2)-1)*blkSize + dim1_2upsample(:), ...
            (index(3)-1)*blkSize + dim2_2upsample(:),...
            points2_out(:, 3:5)]; % x,y,z coord. & color
    elseif idx_var == 2
        blk_return{index(1), index(2), index(3)} = ...
            [(index(1)-1)*blkSize + dim1_2upsample(:),...
            (index(2)-1)*blkSize + dim3_2upsample(:), ...
            (index(3)-1)*blkSize + dim2_2upsample(:),...
            points2_out(:, 3:5)]; % x,y,z coord. & color
    elseif idx_var == 3
        blk_return{index(1), index(2), index(3)} = ...
            [(index(1)-1)*blkSize + dim1_2upsample(:),...
            (index(2)-1)*blkSize + dim2_2upsample(:), ...
            (index(3)-1)*blkSize + dim3_2upsample(:),...
            points2_out(:, 3:5)]; % x,y,z coord. & color
    end
    
end
upsampled_pcloud = [];
% Un-do block encapsulation
for xx = 1 : size(blk_return, 1)
    for yy = 1 : size(blk_return, 2)
        for zz = 1 : size(blk_return, 3)
            upsampled_pcloud = [upsampled_pcloud; blk_return{xx, yy, zz}];
            
        end
    end
end
pcloud_return = pointCloud(upsampled_pcloud(:, 1:3), 'Color', upsampled_pcloud(:, 4:6)/255);
end