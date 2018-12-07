clc;
clear all;
close all;

%% detection
img_path = 'imgs/fig6_b.bmp';
img = imread(img_path);
thresh = graythresh(img);
I2 = im2bw(img, thresh);

T1 = 240;
T2_abs = 190;
T2_rel = 1.2;
N_min = 3000;
T3 = 5;
specular_mask = SpecularDetectionArnold2010(img, T1, T2_abs, T2_rel, N_min, T3);
specular_mask = specular_mask & I2; % filter the background
figure;imshowpair(img,specular_mask,'montage');

%% inpainting
decay_win_size = 10;
decay_cof = 20;
inpainted_img = InpainttingArnold2010(specular_mask, img, decay_win_size, decay_cof);
figure;imshowpair(img,inpainted_img,'montage');
