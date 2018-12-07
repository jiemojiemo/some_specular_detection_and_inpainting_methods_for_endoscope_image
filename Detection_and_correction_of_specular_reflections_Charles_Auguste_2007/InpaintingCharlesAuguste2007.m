function [out] = InpaintingCharlesAuguste2007(img, mask, tol,show_inpainting)
    if nargin <= 2
       tol = 0.05;
       show_inpainting = true;
    end
    
    if nargin <= 3
       show_inpainting = true;
    end

    [nrows ncols nch] = size(img);
    % Checking if the image and the mask have the same size.
    assert(isequal([nrows ncols], size(mask)),"The image and mask must have the same size");
    % Diffusion kernels.
    a = .073235;
    b = .176765;
    c = .125;
    kernel1 = [a b a; b 0 b; a b a];
    kernel2 = [c c c; c 0 c; c c c];
    % Eventual multiple dimensional mask to help some operations with the image's
    % roi.
    mmask = repmat(mask, [1 1 nch]);
    mmask = uint8(mmask);
    mmaskIdxs = find(mmask);
    % Removing color information of the mask pixels from the image.
    out = img .* (1 - mmask);
    % Diffusion iteration.
    maskDiff = Inf;
    
    if show_inpainting
            figure;
    end
    
    while maskDiff > tol
        tempImg = imfilter(out, kernel1, "conv");
        
        % Calculating the roi's difference between 2 successive iterations.
        % We need to use something bigger than the int8 of the images
        % because we need to handle possible negative numbers here.
        maskDiff = mean(abs(int16(out(mmaskIdxs)) - int16(tempImg(mmaskIdxs))));
        out(mmaskIdxs) = tempImg(mmaskIdxs);
        
        if show_inpainting
            imshow(out);
        end
        
    end
end