function [ err ] = gs_errfunction( P0, Xobs )
%GS_ERRFUNCTION: This function calculate the error that we need to pass to the lsqnonlin function
% NOTE: gs_errfunction should return E(X) and not the sum-of-squares E=sum(E(X).^2)) that we want to minimize. 
% (E(X) is summed and squared implicitly in the lsqnonlin algorithm.)
%   Arguments:
%       - P0: Contains the first 9 rows of the homography and then X
%       - Xobs: contains x and x'

% First get the H matrix
    H = reshape(P0(1:9), [3,3]);

    % get x and x' contained in Xobs
    n_points = size(Xobs,1) / 2;
    x = Xobs(1:n_points); % x
    x = reshape(x, [2,size(x,1)/2]);
    xp = Xobs(n_points+1:end); % x'
    xp = reshape(xp, [2,size(xp,1)/2]);
    
    %Error
    xhat = P0(9+1:end);
    xhat = reshape(xhat, [2,size(xhat,1)/2]);
    xhat = [xhat ; ones(1,size(xhat,2))]; %Homogeneous coordinates
    xhatp = H*xhat;
    err = l2_dist(x-euclid(xhat))+l2_dist(xp-euclid(xhatp));

end

function distance = l2_dist(vec)
    distance = (sum(vec.^2));
end

