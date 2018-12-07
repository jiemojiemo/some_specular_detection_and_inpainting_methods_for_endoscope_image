%% HARMONIC Inpainting 
% part of "MATLAB Codes for the Image Inpainting Problem"
%
% Usage:
% u = inpainting_harmonic(u,mask,lambda,tol,maxiter,dt)
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

function u = inpainting_harmonic(input,mask,lambda,tol,maxiter,dt)

[M,N,C] = size(input);

%% LAPLACIAN
% GRID INTERVAL FOR AXIS ij
h1 = 1; h2 = 1;

d2i = toeplitz(sparse([1,1],[1,2],[-2,1]/h1^2,1,M));
d2j = toeplitz(sparse([1,1],[1,2],[-2,1]/h2^2,1,N));
% NEUMANN BOUNDARY CONDITIONS
d2i(1,[1 2])         = [-1 1]/h1;
d2i(end,[end-1 end]) = [1 -1]/h1;
d2j(1,[1 2])         = [-1 1]/h2;
d2j(end,[end-1 end]) = [1 -1]/h2;
% 2D domain LAPLACIAN
L = kron(speye(N),d2i)+kron(d2j,speye(M));

%% ------------------------------------------------------------ FREE MEMORY
clear d2i d2j

%% INPAINTINH ALGORITHM

% INITIALIZATION
u = input;
f = input;

% FOR EACH COLOR CHANNEL
for c = 1:C
    
    for iter = 1:maxiter
        
        % COMPUTE NEW SOLUTION
        laplacian = reshape( L*reshape(u(:,:,c),[],1) ,M,N);
        unew      = u(:,:,c) + dt*( laplacian + lambda*mask(:,:,c).*(f(:,:,c)-u(:,:,c)) );
        
        % COMPUTE EXIT CONDITION
        diff = norm(unew(:)-reshape(u(:,:,c),[],1))/norm(unew(:));
        
        % UPDATE
        u(:,:,c) = unew;
        
        % TEST EXIT CONDITION
        if diff<tol
            break
        end
        
    end  
end

%% WRITE IMAGE OUTPUT
imwrite(u,'./results/harmonic_output.png')

return