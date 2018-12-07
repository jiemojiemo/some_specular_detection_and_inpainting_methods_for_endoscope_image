function [specular_mask] = SpecularDetectionArnold2010(img, T1, T2_abs, T2_rel, N_min, T3)
%SPECULARDETECTIONARNOLD2010 此处显示有关此函数的摘要
%   此处显示详细说明
    cR = double(img(:,:,1)); 
    cG = double(img(:,:,2));
    cB = double(img(:,:,3));
    cE = 0.2989*cR + 0.5870*cG + 0.1140*cB;

    % module 1
    module1_specular_mask = calc_module1_specular_mask(cE, cG, cB, T1);
    
    % module 2
    specular_mask_T2_abs = calc_module1_specular_mask(cE, cG, cB, T2_abs);
    filled_img = filling_image_using_centroid_color(specular_mask_T2_abs, img);
    module2_specular_mask = calc_modul2_specular_mask(filled_img, T2_rel, cR, cG, cB);
    
    final_mask = module1_specular_mask | module2_specular_mask;
    N_min = 5000;
    T3 = 5;
    specular_mask = postprocessing(final_mask, cE, N_min, T3);
end

function [module1_specular_mask] = calc_module1_specular_mask(cE, cG, cB, T1)
    p95_cG = prctile(cG(:), 95);
    p95_cB = prctile(cB(:), 95);
    p95_cE = prctile(cE(:), 95);

    rGE = p95_cG/p95_cE;
    rBE = p95_cB/p95_cE;
    module1_specular_mask = (cG > rGE*T1) | (cB > rBE*T1) | (cE > T1);
end

function module2_specular_mask = calc_modul2_specular_mask(filled_img, T2_rel, cR, cG, cB)
    fR = double(medfilt2(filled_img(:,:,1), [30 30], 'symmetric'));
    fG = double(medfilt2(filled_img(:,:,2), [30 30], 'symmetric'));
    fB = double(medfilt2(filled_img(:,:,3), [30 30], 'symmetric'));
    filtered_img = cat(3, fR, fG, fB);
    
    fR(fR < eps) = 1e7;
    fG(fG < eps) = 1e7;
    fB(fB < eps) = 1e7;
    
    % contrast coefficient
    tR = contrast_coeffcient(single(cR));%tR = 1;
    tG = contrast_coeffcient(single(cG));%tG = 1;
    tB = contrast_coeffcient(single(cB));%tB = 1;
    
    max_img = cat(3, tR*cR./fR, tG*cG./fG, tB*cB./fB);
    e_max = max(max_img, [], 3);
    module2_specular_mask = e_max > T2_rel;
end

function [post_specular_mask] = postprocessing(final_specular_mask, cE, N_min, T3)

    %grad_img = imgradient(cE);
    final_specular_mask = imdilate(final_specular_mask, strel('diamond', 1));
    labeled_area = bwlabel(final_specular_mask);
    num_region = max(labeled_area(:));
    post_specular_mask = final_specular_mask;

    for i=1:num_region
        index = find(labeled_area == i);

        if(length(index) >= N_min)
            post_specular_mask(index) = 0;
        end
    end
end

function [t] = contrast_coeffcient(c)
    mean_c = mean(c);
    std_c = std(c);
    t = 1/((mean_c + std_c)/mean_c);
end