clc;
clear all;
close all;
dbstop if error

addpath ./lib

%% Detection
img_path = 'imgs/2.bmp';
%img_path = 'imgs/image-0001615.png';
img = imread(img_path);
specular_mask = SpecularDetectionOthmane2011(img);
%dilated_mask = imdilate(specular_mask, strel("diamond", 1));
figure;imshow([specular_mask]);

%% Inpaiting
result = inpainting(img, specular_mask, 0.05);
