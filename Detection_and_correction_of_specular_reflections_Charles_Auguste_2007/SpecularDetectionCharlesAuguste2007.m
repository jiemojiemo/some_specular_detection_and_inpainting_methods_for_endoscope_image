function [specular_mask] = SpecularDetectionCharlesAuguste2007(img,inputArg2)

% Step 1: reflection enhancemnt
enhanced = reflection_enhance(im2double(img));
enhanced = im2uint8(enhanced);
enhanced_gray = rgb2gray(enhanced);
figure;imshow([img enhanced]);

% Step 2: histogram denoising
denoised_hist = histogram_denoise(enhanced_gray);

% Step 3: Specular bump thresholding
threshold = find_specular_bump_threshold(denoised_hist);
threshold = 190;

% Step 4: Specular lobe detection
specular_mask_rgb = (enhanced_gray >= (threshold));
%specular_spike_mask = specular_mask_rgb(:,:,1) | specular_mask_rgb(:,:,2) | specular_mask_rgb(:,:,3);
specular_spike_mask = specular_mask_rgb;

specular_mask = imdilate(specular_spike_mask, strel("diamond", 1));

end

function [enhanced] = reflection_enhance(img)
    hsv = rgb2hsv(img);
    s = hsv(:,:,2);
    enhanced = (1-s).*img;
end

function [denoised_hist] = histogram_denoise(enhanced)
    % enhanced is a RGB image and it is double type
    
    [p, x] = imhist(enhanced);
    
    lev = 8;
    denoised_hist = wdenoise(p, lev, 'DenoisingMethod', 'UniversalThreshold');
    
%     figure;
%     plot(x, p);hold on;
%     plot(x, denoised_hist);
end

function [threshold] = find_specular_bump_threshold(w)
    
    w2 = gradient(w);                             
    w2 = thresholding(w2);                            % gradient of the histgram                 
    
    w3 = gradient(w2);
    w3 = thresholding(w3); 
    threshold = find(w3 > 0, 1, 'last');              % threshod is the last positive value
    
end

function [w] = thresholding(w)
    w( w <= 0 ) = 0;
    w( w > 0 ) = 1;
end
