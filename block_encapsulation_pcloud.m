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

% Function that encapsulates the point cloud into dices

function [blks_all, processingOrder] = ...
    block_encapsulation_pcloud(pcloud, blkSize, borderWidth, rho)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:   blks_all    - structure containing all blocks as struct
%                       - each block contains [x, y, z, val, weight]
%           processingOrder - blocks containing more samples go first

% Input:    pcloud - point cloud with locations and color
%           blkSize - for block-wise processing
%           borderWidth - width of the support area in px for point clouds
%                       set to 0
%           rho - decaying factor for weighting function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = double(pcloud.Location(:, 1));
y = double(pcloud.Location(:, 2));
z = double(pcloud.Location(:, 3));
v = double(pcloud.Color);
if isempty(v)
    v = double(pcloud.Intensity);
end

% Remove Offset
x = x - pcloud.XLimits(1);
y = y - pcloud.YLimits(1);
z = z - pcloud.ZLimits(1);


no_blks_x = ceil(diff(pcloud.XLimits)/blkSize);
no_blks_y = ceil(diff(pcloud.YLimits)/blkSize);
no_blks_z = ceil(diff(pcloud.ZLimits)/blkSize);
blks_single = cell(no_blks_x, no_blks_y, no_blks_z);
blks_all = cell(no_blks_x, no_blks_y, no_blks_z);
score = zeros(no_blks_x, no_blks_y, no_blks_z);    % Sets the processing order

% Sort points into blocks
for ii = 1 : length(x) % loop through points
    xblock = ceil(x(ii) / blkSize);
    yblock = ceil(y(ii) / blkSize);
    zblock = ceil(z(ii) / blkSize);
    if xblock == 0
        xblock = 1;
    end
    if yblock == 0
        yblock = 1;
    end
    if zblock == 0
        zblock = 1;
    end
    if xblock <= 0 || yblock <= 0 || zblock<=0 || xblock > no_blks_x || yblock > no_blks_y || zblock > no_blks_z
        continue
    end
    blks_single{xblock, yblock, zblock} = [blks_single{xblock, yblock, zblock}; [x(ii), y(ii), z(ii), v(ii, :)]];
end


disp('--- Point sorting into blocks finished! ---')


% Add support area to points
for yblock = 1 : no_blks_y
    for xblock = 1 : no_blks_x
        for zblock = 1 : no_blks_z
            
            no_neigh_blocks = round(borderWidth / blkSize);
            
            infLimY = max(1, yblock - no_neigh_blocks);
            supLimY = min(no_blks_y, yblock + no_neigh_blocks);
            infLimX = max(1, xblock - no_neigh_blocks);
            supLimX = min(no_blks_x, xblock + no_neigh_blocks);
            infLimZ = max(1, zblock - no_neigh_blocks);
            supLimZ = min(no_blks_z, zblock + no_neigh_blocks);
            
            
            for jj = infLimX : supLimX % x
                for ii = infLimY : supLimY % y
                    for kk = infLimZ : supLimZ % z
                        if isempty(blks_single{jj, ii, kk})
                            continue
                        end
                        xx = blks_single{jj, ii, kk}(:,1) - (infLimX-1) * blkSize;
                        yy = blks_single{jj, ii, kk}(:,2) - (infLimY-1) * blkSize;
                        zz = blks_single{jj, ii, kk}(:, 3) - (infLimZ - 1) * blkSize;
                        blks_all{xblock, yblock, zblock} = ...
                            [blks_all{xblock, yblock, zblock}; [xx, yy, zz, blks_single{jj, ii, kk}(:, 4:end)]];
                    end
                end
            end
            
            if isempty(blks_all{xblock, yblock, zblock})
                score(xblock, yblock, zblock) = -1;
                continue
            end
            
            
            
            % Compute weighting function
            M = (supLimY - infLimY + 1) * blkSize; % y
            N = (supLimX - infLimX + 1) * blkSize; % x
            O = (supLimZ - infLimZ + 1) * blkSize; % z
            xx = blks_all{xblock, yblock, zblock}(:, 1);
            yy = blks_all{xblock, yblock, zblock}(:, 2);
            zz = blks_all{xblock, yblock, zblock}(:, 3);
            val = blks_all{xblock, yblock, zblock}(:, 4:end);
            w = ones(size(xx));
            
            
            for ii = 1:length(xx)
                
                if M < (2*borderWidth + blkSize) && N < (2*borderWidth + blkSize) && O < (2*borderWidth + blkSize)
                    w(ii) = w(ii) * rho^(sqrt(...
                        (yy(ii)-1+0.5-blkSize/2)^2 + ...
                        (xx(ii)-1+0.5-blkSize/2)^2+ ...
                        (zz(ii)-1+0.5-blkSize/2)^2 ...
                        ));
                elseif M < (2*borderWidth + blkSize) && N < (2*borderWidth + blkSize)
                    w(ii) = w(ii) * rho^(sqrt(...
                        (yy(ii)-1+0.5-blkSize/2)^2 + ...
                        (xx(ii)-1+0.5-blkSize/2)^2 + ...
                        (zz(ii)-1+0.5-O/2)^2 ...
                        ));
                    
                elseif M < (2*borderWidth + blkSize) && O < (2*borderWidth + blkSize)
                    w(ii) = w(ii) * rho^(sqrt(...
                        (yy(ii)-1+0.5-blkSize/2)^2 + ...
                        (xx(ii)-1+0.5-N/2)^2 + ...
                        (zz(ii)-1+0.5-blkSize/2)^2 ...
                        ));
                    
                elseif N < (2*borderWidth + blkSize) && O < (2*borderWidth + blkSize)
                    w(ii) = w(ii) * rho^(sqrt(...
                        (yy(ii)-1+0.5-M/2)^2 + ...
                        (xx(ii)-1+0.5-blkSize/2)^2 + ...
                        (zz(ii)-1+0.5-blkSize/2)^2 ...
                        ));
                elseif M < (2*borderWidth + blkSize)
                    w(ii) = w(ii) * rho^(sqrt(...
                        (yy(ii)-1+0.5-blkSize/2)^2 + ...
                        (xx(ii)-1+0.5-N/2)^2 + ...
                        (zz(ii)-1+0.5-O/2)^2 ...
                        ));
                elseif N < (2*borderWidth + blkSize)
                    w(ii) = w(ii) * rho^(sqrt(...
                        (yy(ii)-1+0.5-M/2)^2 + ...
                        (xx(ii)-1+0.5-blkSize/2)^2 + ...
                        (zz(ii)-1+0.5-O/2)^2 ...
                        ));
                elseif O < (2*borderWidth + blkSize)
                    w(ii) = w(ii) * rho^(sqrt(...
                        (yy(ii)-1+0.5-M/2)^2 + ...
                        (xx(ii)-1+0.5-N/2)^2 + ...
                        (zz(ii)-1+0.5-blkSize/2)^2 ...
                        ));
                else
                    w(ii) = w(ii) * rho^(sqrt(...
                        (yy(ii)-1+0.5-M/2)^2 + ...
                        (xx(ii)-1+0.5-N/2)^2 + ...
                        (zz(ii)-1+0.5-O/2)^2 ...
                        ));
                end
                
                
            end
            blks_all{xblock, yblock, zblock} = [blks_all{xblock, yblock,zblock}, w];
            
            
            % Processing order
            if all(val(:,1) < 0) || isempty(val) % R2017b
                score(xblock, yblock, zblock) = -1;
            else
                score(xblock, yblock, zblock) = sum(w);
            end
            
            
        end
    end
end

disp('--- Support Area finished! ---')

% Obtaining processing order
[yidx, xidx, zidx] = meshgrid(1:no_blks_y, 1:no_blks_x, 1:no_blks_z);
score = score(:); % mapping gets lost
yidx = yidx(:);
xidx = xidx(:);
zidx = zidx(:);
[~, idx] = sort(score, 'descend');      % The more samples inside the better
numBlocks = sum(score >= 0);            % Empty blocks are not processed
yidx = yidx(idx);
xidx = xidx(idx);
zidx = zidx(idx);
processingOrder = [xidx(1:numBlocks)'; yidx(1:numBlocks)'; zidx(1:numBlocks)'];

disp('--- Processing Order finished! ---')

end



