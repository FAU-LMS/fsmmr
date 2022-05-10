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

function [points2_in, points2_out] = ...
    transformSurface2Plane_shortestWay(points3_in, points3_out, trafo_size, RGB_tmp)

%% Transform for the final/ resampled  points
points3_all = [points3_in; points3_out];
dim1 = points3_all(:,1);
dim2 = points3_all(:,2);
dim3 = points3_all(:,3);
dim1_mat = dim1 - dim1';
dim2_mat = dim2 - dim2';
dim3_mat = dim3 - dim3';
% dim1_mat(dim1_mat==0) = 1000;
% dim2_mat(dim2_mat==0) = 1000;
% dim3_mat(dim3_mat==0) = 1000;
euklid_mat = sqrt(dim1_mat.^2+dim2_mat.^2+dim3_mat.^2);
[min_euklid, idx] = min(euklid_mat);

G = graph(euklid_mat);
[T, PRED] = minspantree(G);

delta1_mat = sign(dim1_mat) .* double(sqrt(dim1_mat.^2+dim3_mat.^2));
delta2_mat = sign(dim2_mat) .* double(sqrt(dim2_mat.^2+dim3_mat.^2));

coord2d = zeros(size(euklid_mat,1), 2);
delta1_start = sqrt(dim1(1).^2 + dim3(1).^2);
delta2_start = sqrt(dim1(2).^2 + dim3(1).^2);
for ii = 2 : size(euklid_mat, 1)
    shortPath = shortestpath(T, 1, ii);
    pathLength1 = 0;
    pathLength2 = 0;
    for ll = 1 : length(shortPath)-1
        pathLength1 = pathLength1 + delta1_mat(shortPath(ll), shortPath(ll+1));
        pathLength2 = pathLength2 + delta2_mat(shortPath(ll), shortPath(ll+1));
    end
    coord2d(ii, 1) = pathLength1;% + delta1_start;
    coord2d(ii, 2) = pathLength2;% + delta2_start;
end

coord2d(:, 1) = coord2d(:, 1) - min(coord2d(:, 1));
coord2d(:, 2) = coord2d(:, 2) - min(coord2d(:, 2));

% figure
% for kk = 1 : length(dim1)
%     if kk <= length(points3_in)
%         plot(coord2d(kk, 1), coord2d(kk, 2), 'x', 'Color', RGB_tmp(kk, :)/255)
%     else
%         plot(coord2d(kk, 1), coord2d(kk, 2), 'o', 'Color', 'b')
%     end
%     hold on
% end

points2_out = [coord2d(size(points3_in, 1)+1 : end, 1), coord2d(size(points3_in, 1)+1 : end, 2)];
points2_in = [coord2d(1:size(points3_in, 1), 1), coord2d(1:size(points3_in, 1), 2)];
end




