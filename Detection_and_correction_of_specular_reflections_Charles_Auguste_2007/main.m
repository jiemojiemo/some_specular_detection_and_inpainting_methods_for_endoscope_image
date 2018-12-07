clc;
clear all;
close all;

dbstop if error;

show_mask = true;

% read image
%img_path = 'imgs\image-0001615.png';
img_path = 'imgs\figure9_2.bmp';
img = imread(img_path);

%% Detection

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

decrease_threshod = calculate_decrease_threshold(img);
win_size = 3;
enlarged_mask = enlarge_spike_mask(img, specular_spike_mask, win_size, decrease_threshod);

dilated_mask = imdilate(specular_spike_mask, strel("diamond", 1));

if show_mask
    figure;
    imshow([specular_spike_mask enlarged_mask dilated_mask]);
end

%% correction
inpaiting_img = inpainting(img, dilated_mask, 0.005);
figure;imshow([img inpaiting_img]);

