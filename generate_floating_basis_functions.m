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

function [basis_functions, frequency_weighting] = generate_floating_basis_functions(y, x, ...
    transform_size, transform_type, frequency_weighting_mode, sigma)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output: 
% basis_functions: dct basis functions evaluated at all x,y positions
% frequency_weighting: frequency weighting determined for all frequencies
% 
% Input: 
% y, x: coordinates in dim1, dim2
% transform_size: chosen transform size for dct transform
% transform_type: here, only 'dct' possible
% frequency_weighting_mode: normally set to '5', '1' corresponds to no
% frequency weighting
% sigma: decaying factor for frequency weighting function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

len = length(y);
basis_functions = zeros(len, transform_size^2);
frequency_weighting_factor = 2;

if strcmp(transform_type, 'dct')
    % discrete cosine transform
    len = length(y);
    % weighting of dct basis functions
    w = sqrt(2/transform_size)*ones(1, transform_size);
    w(1) = sqrt(1/transform_size);
    ww = w'*w;
    
    basis_functions = zeros(len, transform_size^2);
    % iterate through positions
    for n = 1:len
        bf_counter = 1;
        % iterate frequencies
        for k1 = 1:transform_size
            for k2 = 1:transform_size
                basis_functions(n, bf_counter) = ww(k1,k2) * ...
                    cos( pi/(2*transform_size) * (2*(y(n)-1) + 1)*(k1-1) ) ...
                    * cos( pi/(2*transform_size) * (2*(x(n)-1) + 1)*(k2-1) );
                bf_counter = bf_counter + 1;
            end % for k2
        end % for k1
        
    end % for n
    
% frequency weighting
frequency_weighting = ones(transform_size^2, 1);
if(frequency_weighting_mode == 1)
    frequency_weighting(:) = 1;
    
elseif (frequency_weighting_mode == 2)
    bf_counter = 1;
    for y_counter = 0:transform_size-1
        for x_counter = 0:transform_size-1
            y2 = transform_size/2 - abs(y_counter - transform_size/2);
            x2 = transform_size/2 - abs(x_counter - transform_size/2);
            frequency_weighting(bf_counter) = 1 - sqrt(x2*x2 + y2*y2)*sqrt(2)/transform_size;
            bf_counter = bf_counter + 1;
        end
    end
    
elseif (frequency_weighting_mode == 3) % bath tub
    bf_counter = 1;
    for y_counter=0:transform_size-1
        for x_counter=0:transform_size-1
            y2 = transform_size/2 - abs(y_counter - transform_size/2);
            x2 = transform_size/2 - abs(x_counter - transform_size/2);
            frequency_weighting(bf_counter) = exp(-frequency_weighting_factor*sqrt(x2*x2 + y2*y2)*sqrt(2)/transform_size);
            bf_counter = bf_counter + 1;
        end
    end
    
elseif (frequency_weighting_mode == 4)
    for y_counter=0:transform_size-1
        for x_counter=0:transform_size-1
            y2 = transform_size/2 - abs(y_counter - transform_size/2);
            x2 = transform_size/2 - abs(x_counter - transform_size/2);
            frequency_weighting(y_counter+1, x_counter+1) = exp(-frequency_weighting_factor*(sqrt(x2*x2 + y2*y2)*sqrt(2)/transform_size)^2);
        end
    end
        
elseif  (frequency_weighting_mode == 5) % FSMMR
    bf_counter = 1;
    fw = zeros(transform_size);
    for y_counter=0:transform_size-1
        for x_counter=0:transform_size-1
            y2 = y_counter;
            x2 = x_counter;
            frequency_weighting(bf_counter) = sigma^(sqrt(x2*x2 + y2*y2));
            bf_counter = bf_counter + 1;
            fw(y_counter+1,x_counter+1) = frequency_weighting(bf_counter-1);
        end %for x_counter
    end %for y_counter    
end %if frequency_weighting_mode

end %function