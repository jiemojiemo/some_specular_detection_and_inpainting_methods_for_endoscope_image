clc;
clear all;
close all;

%% Detection
img_path = 'imgs/figure8_1.bmp';
img = imread(img_path);
specular_mask = SpecularDetectionCharlesAuguste2007(img);
figure;imshow(specular_mask);

%% Inpaiting
inpaited_img = InpaintingCharlesAuguste2007(img, specular_mask, 0.05);