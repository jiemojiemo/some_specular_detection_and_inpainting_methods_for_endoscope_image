clc;
clear all;
close all;
dbstop if error

addpath ./lib

%% Detection
img_path = 'imgs/fig10_b.bmp';
%img_path = 'imgs/image-0001615.png';
img = imread(img_path);
specular_mask = SpecularDetectionMeslouhi2011(img);
%dilated_mask = imdilate(specular_mask, strel("diamond", 1));
imshowpair(img, specular_mask, 'montage');

%% Inpaiting
%result = inpainting(img, specular_mask, 0.05);
