function [inpainted_img] = InpainttingArnold2010(specular_mask,img, decay_win_size, decay_cof)
    filled_img = filling_image_using_centroid_color(specular_mask, img);
    % gaussian filter
    sig = 8;
    gaussian_filtered_img = imgaussfilt(filled_img, sig);

    % mask decay (not the paper solution but very similar to it)
    mx = imfilter(double(specular_mask), ones(decay_win_size)/decay_cof);
    mx = mx + specular_mask;
    mx(mx > 1) = 1.0;
    inpainted_img = mx.*double(gaussian_filtered_img) + (1-mx).*double(img);
    inpainted_img = uint8(inpainted_img);
end

