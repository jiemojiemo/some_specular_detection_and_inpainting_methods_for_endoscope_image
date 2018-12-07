function [filled_img] = filling_image_using_centroid_color(specular_mask_T2_abs, img)

    centroid_color_info = calc_centroid_color_info(specular_mask_T2_abs, img);
    
    % filling possible specular highlights by the centroid color
    [row_index, col_index] = find(specular_mask_T2_abs == 1);
    num_possible_specular_points = length(row_index);
    filled_img = img;
    for i=1:num_possible_specular_points
        % looking for the nearst centroid color for every specular point
        % and fill it
        nearest_region = find_the_nearest_region(centroid_color_info, row_index(i), col_index(i));
        filled_img(row_index(i), col_index(i), :) = nearest_region.centroid_color;
    end
 end
function [centroid_color_info] = calc_centroid_color_info(specular_mask_T2_abs, img)
    
    dilated_mask_1 = imdilate(specular_mask_T2_abs, strel('disk', 2));
    dilated_mask_2 = imdilate(specular_mask_T2_abs, strel('disk', 4));
    centroid_color_area = dilated_mask_2-dilated_mask_1;
    labeled_area = bwlabel(centroid_color_area);
    num_region = max(labeled_area(:));
    
    for i=1:num_region
        [row_index, col_index] = find(labeled_area == i);
        centroid_color_info(i).centroid_row = mean(row_index);
        centroid_color_info(i).centroid_col = mean(col_index);
        centroid_color_info(i).centroid_color = mean(img(row_index, col_index, :), [1 2]);
    end
end
function nearest_region = find_the_nearest_region(centroid_color_info, pixel_row, pixel_col)

    num_region = length(centroid_color_info);
    nearset_region_index = 0;
    nearset_distance = 1e6;
    for j=1:num_region
            distance_to_centroid = calc_distance( pixel_row, pixel_col,...
                centroid_color_info(j).centroid_row, centroid_color_info(j).centroid_col);
            
            if(distance_to_centroid < nearset_distance)
                nearset_distance = distance_to_centroid;
                nearset_region_index = j;
            end
    end
    
    nearest_region = centroid_color_info(nearset_region_index);
end
function [distance_to_centroid] = calc_distance(x,y,x1,y1)
    distance_to_centroid = sqrt( (x-x1)^2 + (y-y1)^2 );
end