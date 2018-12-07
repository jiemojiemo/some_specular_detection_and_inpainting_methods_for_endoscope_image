function [enlarged_mask] = enlarge_spike_mask(img, specular_spike_mask, win_size, threshold)
    decrease_threshod = calculate_decrease_threshold(img);
    
    [num_rows, num_cols] = size(specular_spike_mask);
    
    [spike_row, spike_col] = find(specular_spike_mask);
    
    num_spike = length(spike_row);
    
    enlarged_mask = specular_spike_mask;
    for i=1:num_spike
       center_row = spike_row(i);
       center_col = spike_col(i);
       
       win_location_i_begin = center_row - floor(win_size/2);
       if(win_location_i_begin < 1) win_location_i_begin = 1; end
       
       win_location_i_end = center_row + floor(win_size/2);
       if(win_location_i_end > num_rows) win_location_i_end = num_rows; end
       
       win_location_j_begin = center_col - floor(win_size/2);
       if(win_location_j_begin < 1) win_location_j_begin = 1; end
       
       win_location_j_end = center_col + floor(win_size/2);
       if(win_location_j_end > num_cols) win_location_j_end = num_cols; end
       
       for j=win_location_i_begin:win_location_i_end
          for k=win_location_j_begin:win_location_j_end
              distance = sqrt( (j-center_row)^2 + (k-center_col)^2 );
              new_threshold = threshold - distance*decrease_threshod;
              
              temp = img(center_row, center_col, :) > new_threshold;
              is_specular_lobe = temp(:,:,1) | temp(:,:,2) | temp(:,:,3);
              
              if(is_specular_lobe)
                  enlarged_mask(j, k) = 1;
              end
          end
       end
       
    end
    
    %[num_row, num_col] = size(specular_spike_mask);
end