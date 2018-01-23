function [disparity, matching_cost] = stereo_computation(left_img, right_img, min_disp, max_disp, window_size, cost)

    
    disparity = zeros(size(left_img));    
    [h,w] = size(left_img);
    %Addition of padding for sliding window.
    padding = floor(window_size / 2);
    left_img = padarray(left_img, [padding padding]);
    right_img = padarray(right_img, [padding padding]);
    
    for x=1+padding:h+padding
        for y=1+padding:w+padding
            left = left_img(x-padding:x+padding, y-padding:y+padding);
            min_left = max(1+padding, y-max_disp);
            max_left = min(w+padding, y+max_disp);
            weights = zeros(size(left));
            weights(:) = 1/numel(left);
            if cost=='SSD'
                min_ssd = Inf; %Cost function --> minimize
                for i = min_left:max_left
                   right = right_img(x-padding:x+padding, i-padding:i+padding); 
                   ssd_value = sum(sum(weights.*(left-right).^2 ));
                   if ssd_value < min_ssd
                        min_ssd = ssd_value;
                        index_best_cost = i;
                    end
                end
                matching_cost = min_ssd;
            elseif cost == 'NCC'
                max_ncc = -Inf;  %Similarity metric --> maximize             
                for i = min_left:max_left
                    right = right_img(x-padding:x+padding, i-padding:i+padding);
                    sum_left = sum(left(:).*weights(:));
                    sum_right = sum(right(:).*weights(:));
                    
                    sigma_left = sqrt(sum(weights(:).* (left(:) - sum_left).^2 ));
                    sigma_right = sqrt(sum(weights(:).* (right(:) - sum_right).^2 ));
                    
                    ncc_value = sum(weights(:).*(left(:)-sum_left).*(right(:)-sum_right))/(sigma_left*sigma_right);
                    if ncc_value > max_ncc
                        max_ncc = ncc_value;
                        index_best_cost = i;
                    end
                end
                matching_cost = max_ncc;
            end
            disparity(x-padding, y-padding) = abs(y-index_best_cost);
        end
    end


end