function [Pproj, Xproj] = factorization_method(x1, x2, sturm)

    if nargin < 3
        sturm = false;
    end
    d = inf;
    num_points = size(x1,2)
    %Normalize set of points in image:
    [x1_norm, T1] = normalise2dpts(x1);
    [x2_norm, T2] = normalise2dpts(x2);

    %Initialize lambda:
    if ~sturm
         lambda_matrix = ones (2, num_points);
    else
        %Compute epipolar lines:
        F1 = fundamental_matrix(x1, x1);
        [U, D, V] = svd(F1);
        epipolar_1 = V(:,3) / V(3,3);
        F2 = fundamental_matrix(x2, x1);
        [U, D, V] = svd(F2);
        epipolar_2 = V(:,3) / V(3,3);
        
        %Allocate memory for lambda matrix
        lambda_matrix = ones (2, num_points);
        for i=1:size(x1,2)
           lambda_matrix(1, i) =  x1(:, i)'*F1*cross(epipolar_1, x1(:,i)) ...
               / norm(cross(epipolar_1, x1(:,i))).^2;
           lambda_matrix(2, i) =  x1(:, i)'*F2*cross(epipolar_2, x2(:,i)) ...
               / norm(cross(epipolar_2, x2(:,i))).^2;
        end
    end

    while 1
       dist = Inf;
       iters = 0;
       %Alternate rescaling row and columns no have unit norm
       while 1
           previous_dist = dist;
           previous_lambda = lambda_matrix
           if mod(iters,2) > 0
               lambda_matrix(1,:) = lambda_matrix(1,:) ./ norm(lambda_matrix(1,:));
               lambda_matrix(2,:) = lambda_matrix(2,:) ./ norm(lambda_matrix(2,:));
               
           else
                for j = 1:num_points
                    lambda_matrix(:,j) = lambda_matrix(:,j) / norm(lambda_matrix(:,j));
                end
           end
           dist = (previous_lambda-lambda_matrix).^2;
           dist = sqrt(sum(dist(:)));
           if (abs(dist - previous_dist) / dist) < 0.1
              break 
           end
           iters = iters + 1;
       end
       
       %Build Measurement matrix M:
       M = zeros(6, num_points);
       M(1,:) = lambda_matrix(1,:) .* x1_norm(1,:);
       M(2,:) = lambda_matrix(1,:) .* x1_norm(2,:);
       M(3,:) = lambda_matrix(1,:) .* x1_norm(3,:);
       M(4,:) = lambda_matrix(2,:) .* x2_norm(1,:);
       M(5,:) = lambda_matrix(2,:) .* x2_norm(2,:);
       M(6,:) = lambda_matrix(2,:) .* x2_norm(3,:);
        
       %SVD of M:
       [U,D,V] = svd(M);
       Xproj = V(:,1:4)';
       Pmotion = U*D(:,1:4);
       
       previous_d = d; d = 0;
       P_x1 = Pmotion(1:3,:) * Xproj;
       P_x2 = Pmotion(4:6,:) * Xproj;
       for k=1:num_points
           d = d + sum((x1_norm(:,k) - P_x1(:,k)).^2) + ...
               sum((x2_norm(:,k) - P_x2(:,k)).^2);
       end
       if (abs(d - previous_d) / d) < 0.1
          break
       else
           %Not converged yet:
           P = Pmotion * Xproj;
           lambda_matrix(1,:) = P(3,:);
           lambda_matrix(2,:) = P(6,:);
       end
    end
    Pproj(1:3,:) = inv(T1) * Pmotion(1:3,:);
    Pproj(4:6,:) = inv(T2) * Pmotion(4:6,:);
end