function M = CBIG_uniform_rand_rotation(num_rotations, seed)

% M = CBIG_uniform_rand_rotation(num_rotations, seed)
%
% This function is used to compute a uniformly distributed rotation matrix
% (UDRM). To compute the UDRM, and from a macro/big-picture view, the
% function performs the fast random rotation matrix method. This involves
% generating the householder matrix H, and the simple rotation matrix R to
% compute M, the UDRM.
%
% INPUT:
%       - num_rotations
%         A scalar value used to indicate the desired number of rotations.
%         The default value is 1000.
%
%       - seed
%         A scalar value used to define the random seed for the purpose of
%         producing random rotation.
%         The default value is 1.
%
% OUTPUT:
%       - M:
%         A 3x3xnum_rotations matrix, the UDRM for each rotation
%
% Example:
% num_rotations = 5000;
% seed = 2;
% M = CBIG_uniform_rand_rotation(num_rotations, seed);
%Written by Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    if nargin == 0 % If no input for num_rotations, set default num_rotations as 1000
        num_rotations = 1000;
    end

    if nargin <= 1 % If no input for seed, set default seed as 1
        seed = 1;
    end
    
    rng(seed,'twister'); % Setting the rng using the input seed, with Mersenne Twister random number algorithm
    
    for i = 1:num_rotations

        randnums = rand(1,3); % Randomly generating a row vector of 3 double values within range [0 1]: x1, x2, and x3

        theta = 2*pi*randnums(1); % Defining random convenience variable theta based on random x1 value
        phi = 2*pi*randnums(2); % Defining random convenience variable phi based on random x2 value
        x3 = randnums(3); % Defining x3
        r = sqrt(randnums(3)); % Defining random convenience variable r based on square root of random x3 value
        
        R = [cos(theta) sin(theta) 0; % R is a matrix for simple rotation
            -1*sin(theta) cos(theta) 0;
            0 0 1]; % Note that row 3 of R: [0 0 1] represents z, the North Pole
        
        % Defining V which will allow H to satisfy condition of z under the random reflection
        % is such that both its animuthal angle and the cosine of its polar
        % angles are uniformly distributed
        V = [ cos(phi)*r; sin(phi)*r; sqrt(1-x3) ]; 
        
        H = (eye(3) - 2*V*V'); % Defining the Householder matrix
        
        M(:,:,i) = -H*R; % Computing the final uniform rotation matrix
    end
end