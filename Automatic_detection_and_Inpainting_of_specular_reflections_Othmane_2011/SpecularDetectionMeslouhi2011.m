function [specular_mask] = SpecularDetectionMeslouhi2011(img)
    
    % Step 1 Imange Enhancement
    enhanced = ReflectionEnhance(im2double(img));
%     enhanced = im2uint8(enhanced);
%     imshow([img enhanced]);

    % Step 2 Convertion to CIE-XYZ and then getting the lunminace(Y)
    enhanced_xyz = rgb2xyz(enhanced);
    Y = enhanced_xyz(:,:,2);
    
    % Step 3 Getting chromatic luminance(y)
    y = enhanced_xyz(:,:,2)./(enhanced_xyz(:,:,1) + enhanced_xyz(:,:,2) + enhanced_xyz(:,:,3));
    
    % Step 3 Getting Mask
    specular_mask = Y >= y;
    %dilated_mask = imdilate(specular_mask, strel("diamond", 1));
end

function [enhanced] = ReflectionEnhance(img)
    min_rgb = min(img, [], [3]);
    max_rgb = max(img, [], [3]);
    s = min_rgb ./ max_rgb;
    enhanced = s.*img;
end

