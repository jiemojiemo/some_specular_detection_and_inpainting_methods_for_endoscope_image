%% Create data to be inpainted
% part of "MATLAB Codes for the Image Inpainting Problem"
%
% Usage:
% [u,mask,input] = create_image_and_mask(imagefilename,maskfilename)
%
% Authors:
% Simone Parisotto          (email: sp751 at cam dot ac dot uk)
% Carola-Bibiane Schoenlieb (email: cbs31 at cam dot ac dot uk)
%      
% Address:
% Cambridge Image Analysis
% Centre for Mathematical Sciences
% Wilberforce Road
% CB3 0WA, Cambridge, United Kingdom
%  
% Date:
% September, 2016
%
% Licence: BSD-3-Clause (https://opensource.org/licenses/BSD-3-Clause)
%


function [u,mask,input] = create_image_and_mask(imagefilename,maskfilename)

% import a clean input to be corrupted with the mask 
input = im2double(imread(imagefilename));

% import the mask of the inpainting domain
% mask = 1 intact part
% mask = 0 missing domain
mask  = double( mat2gray( im2double(imread(maskfilename)) ) == 1 );
if size(mask,3)==1 && size(input,3)>1
    mask = repmat(mask,[1,1,size(input,3)]);
end

% create the image with the missin domain:
noise = rand(size(input));
u     = mask.*input + (1-mask).*noise;