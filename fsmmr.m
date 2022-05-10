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

function used_freqs_2d = ...
    fsmmr(val_up, w_up, max_iter, ...
    orthogonality_correction, basis_functions_up, freq_weights)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%           val_up:               Color information of block point cloud
%           w_up:                 Spatial weights of upsampled points
%           max_iter:               Maximum number of iterations
%           orthogonality_correction: Gamma from FSR, orthogonality deficiency compensation
%                                   -> stable estimation and reduced interference btw basis functions
%           basis_functions_up:        Basis functions for upsampled points
%           freq_weights:           Frequency weighting function
% Output:
%           used_freqs_2d:      Expansion coefficients of all used
%                               frequencies in this block summed up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialisation Test
used_bfs = zeros(max_iter, 2);

% Source signal
f = val_up;                        % Image signal of block
[~, trafo_length_sq] = size(basis_functions_up);      
trafo_length = sqrt(trafo_length_sq);

% Orthogonality matrix
ortho_matrix_num = ((w_up.*basis_functions_up).')*basis_functions_up; 

% First iteration
res1 = f.*w_up;

% Transform
exp_coeff_num = zeros(trafo_length_sq, 1);
for kk = 1:trafo_length_sq
    exp_coeff_num(kk) = (sum( res1 .* conj(basis_functions_up(:,kk)) ));
end

ortho_matrix_denominator = 1 ./ diag(ortho_matrix_num);

% Further iterations
for iter_counter = 1 : max_iter
    
    % Choose appropriate basis function
    delta_resEnergy = (exp_coeff_num.^2) .* ortho_matrix_denominator .* freq_weights;
    
    bf2select = find(delta_resEnergy == max(delta_resEnergy));
    bf2select = bf2select(1)-1;
    
    % Estimated expansion coefficient
    expansion_coefficient = orthogonality_correction * exp_coeff_num(bf2select+1) / ortho_matrix_num(bf2select+1, bf2select+1);
    
    % Update model
    used_bfs(iter_counter, 1) = bf2select;
    if iter_counter == 1
        used_bfs(iter_counter, 2) = expansion_coefficient;
    else 
        if abs(expansion_coefficient) > abs(used_bfs(iter_counter-1, 2))+50
            used_bfs(iter_counter, 2) = sign(expansion_coefficient) * used_bfs(iter_counter-1, 2);
        else
            used_bfs(iter_counter, 2) = expansion_coefficient;
        end
    end
    % Transform
    exp_coeff_num = exp_coeff_num - expansion_coefficient * ortho_matrix_num(:, bf2select+1);
    
end

% Determine transform coefficients
used_freqs = zeros(trafo_length_sq, 1);
for ii = 0 : (trafo_length_sq - 1)
    [idx, ~] = find((used_bfs(:, 1)) == ii);
    used_freqs(ii+1) = sum(used_bfs(idx, 2));
end
% Take transform coefficients in 2D 
used_freqs_2d = reshape(used_freqs, [trafo_length, trafo_length]);

end
