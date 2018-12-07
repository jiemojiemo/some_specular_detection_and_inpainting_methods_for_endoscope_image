function [threshold_rgb] = calculate_decrease_threshold(img)
    % img is uint8 type

    max_rgb = max(img, [], [1,2,3]);
    min_rgb = min(img, [], [1,2,3]);
    
    threshold_rgb = (max_rgb - min_rgb)/100;
    threshold_rgb = squeeze(threshold_rgb);
end