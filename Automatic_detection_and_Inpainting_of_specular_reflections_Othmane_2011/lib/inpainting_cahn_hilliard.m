%% Cahn-Hilliard Inpainting 
% part of "MATLAB Codes for the Image Inpainting Problem"
%
% Usage:
% u = inpainting_cahn_hilliard(u,mask,maxiter,param)
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

function u = inpainting_cahn_hilliard(input,mask,maxiter,param)
% Cahn-Hilliard inpainting based on the paper
%
% Bertozzi, Andrea L., Selim Esedoglu, and Alan Gillette.
% "Inpainting of binary images using the Cahn-Hilliard equation."
% IEEE Transactions on image processing 16.1 (2007): 285-291.
%
% The modified Cahn-Hilliard equation was discretized based on
% convexity splitting proposed in the same paper and analysed in
%
% C.-B. Schoenlieb, A. Bertozzi, Unconditionally stable schemes for
% higher order inpainting, Communications in Mathematical Sciences,
% Volume 9, Issue 2, pp. 413-457 (2011).
%
% Namely
%
% E_1  = \int_{\Omega} \ep/2 |\nabla u(:,:,c)|^2 + 1/\ep W(u(:,:,c)) dx , W(u(:,:,c)) = u(:,:,c)^2 (1-u(:,:,c))^2
% E_11 = \int_{\Omega} \ep/2 |\nabla u(:,:,c)|^2 + C_1/2 |u(:,:,c)|^2 dx, E_12 =
% \int_{\Omega} - 1/\ep W(u(:,:,c)) + C_1/2 |u(:,:,c)|^2 dx
%
% E_2 = \lambda \int_{\Omega\D} (f-u(:,:,c))^2 dx,
% E_21 = \int_{\Omega\D} C_2/2 |u(:,:,c)|^2 dx, E_22 = \int_{\Omega\D} -\lambda (f-u(:,:,c))^2 + C_2/2|u(:,:,c)|^2 dx

%%  IMPORT THE CLEAN INPUT AND THE MASK
[M,N,C] = size(input);
u       = input;

%% PARAMETERS
% GRID INTERVAL FOR AXIS ij
h1 = 1; h2 = 1;
swap         = round(maxiter/2);
ep           = [param.epsilon(1)*ones(numel(1:swap-1),1); param.epsilon(2)*ones(numel(swap:maxiter),1)];
c1           = 1/param.epsilon(2);
lambda       = param.lambda*mask;
dt           = param.dt;

% Diagonalize the Laplace Operator by: Lu + uL => D QuQ + QuQ D, where 
% Q is nonsingular, the matrix of eigenvectors of L and D is a diagonal matrix.
% We have to compute QuQ. This we can do in a fast way by using the fft-transform:

Lambda1 = spdiags(2*(cos(2*(0:M-1)'*pi/M)-1),0,M,M)/h1^2;
Lambda2 = spdiags(2*(cos(2*(0:N-1)'*pi/N)-1),0,N,N)/h2^2;

Denominator = Lambda1*ones(M,N) + ones(M,N)*Lambda2;

% Now we can write the above equation in much simpler way and compute the
% solution u_hat

for c = 1:C
    
    % Initialization of Fourier transform:
    u_hat      = fft2(u(:,:,c));
    lu0_hat    = fft2(lambda(:,:,c).*u(:,:,c));
    
    for it = 1:maxiter
        
        lu_hat     = fft2(lambda(:,:,c).*u(:,:,c));
        Fprime_hat = fft2(2*(2*u(:,:,c).^3-3*u(:,:,c).^2+u(:,:,c)));
        
        % CH-inpainting
        u_hat = (dt*(1+lambda(:,:,c)-Denominator/param.epsilon(2)).*u_hat...
            + dt/ep(it)*Denominator.*Fprime_hat...
            + dt*(lu0_hat-lu_hat))./(1+lambda(:,:,c)*dt+ep(it)*dt*Denominator.^2-dt*Denominator/param.epsilon(2));
        
        u(:,:,c) = real(ifft2(u_hat));
        
    end

end

%% WRITE IMAGE OUTPUTS
imwrite(u,'./results/cahn_hilliard_output.png')
