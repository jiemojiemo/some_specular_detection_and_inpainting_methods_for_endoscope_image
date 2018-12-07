%% Mumford_Shah Inpainting 
% part of "MATLAB Codes for the Image Inpainting Problem"
%
% Usage:
% [u, chi] = inpainting_mumford_shah(u,mask,maxiter,tol,param)
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

function [u, chi] = inpainting_mumford_shah(input,mask,maxiter,tol,param)
% Inpainting with the Mumford-Shah image model and Ambrosio-Tortorelli.
% For a given grey value image ustart with image domain \Omega and
% inpainting  domain (damaged part) D we want to reconstruct an image u from f
% by solving
% %%%%%%%%%%%%%% MINIMISATION PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% u = argmin_u  (\frac{\alpha}{2} \int_\Omega \chi^2 |\nabla u|^2 dx   %
%               + \beta \int_\Omega \left(\epsilon |\nabla \chi|^2     %
%               + \frac{(1-\chi)^2}{4\epsilon} \right) dx              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The above minimisation problem is solved iteratively via alternating
% solutions of the Euler-Lagrange equations for u and \chi.

%% IMPORT THE CLEAN INPUT AND THE MASK
[M,N,C] = size(input);

%% GRADIENT
% GRID INTERVAL FOR AXIS ij
h1 = 1; h2 = 1;
% GRID INTERVAL FOR AXIS ij
%h1 = 1/(size(input,1)+1); h2 = 1/(size(input,2)+1);

% FORWARD AND BACKWARD
d1i_forward  = spdiags([-ones(M,1),ones(M,1)],[0,1],M,M)/h1;
d1j_forward  = spdiags([-ones(N,1),ones(N,1)],[0,1],N,N)/h2;
d1i_backward = spdiags([-ones(M,1),ones(M,1)],[-1,0],M,M)/h1;
d1j_backward = spdiags([-ones(N,1),ones(N,1)],[-1,0],N,N)/h2;
% PERIODIC BOUNDARY CONDITIONS
d1i_forward(end,[1 end]) = [1 -1]/h1;
d1j_forward(end,[1 end]) = [1 -1]/h2;
d1i_backward(1,[1 end]) = [-1 1]/h1;
d1j_backward(1,[1 end]) = [-1 1]/h2;

matrices.Dif  = kron(speye(N),d1i_forward);
matrices.Dib  = kron(speye(N),d1i_backward);
matrices.Djf  = kron(d1j_forward,speye(M));
matrices.Djb  = kron(d1j_backward,speye(M));

% CENTRAL
matrices.Dic = (matrices.Dif+matrices.Dib)/2;
matrices.Djc = (matrices.Djf+matrices.Djb)/2;

%% -------------------------------------------------------------- LAPLACIAN
d2i = toeplitz(sparse([1,1],[1,2],[-2,1],1,M))/h1^2;
d2j = toeplitz(sparse([1,1],[1,2],[-2,1],1,N))/h2^2;
% PERIODIC BOUNDARY CONDITIONS
d2i(1,2)       = 2/h1^2;
d2i(end,end-1) = 2/h1^2;
d2j(1,2)       = 2/h2^2;
d2j(end,end-1) = 2/h2^2;
% 2D domain LAPLACIAN
matrices.LAP = (kron(speye(N),d2i)+kron(d2j,speye(M)));

%% ------------------------------------------------------------ FREE MEMORY
clear d1i_forward dji_forward d1i_backward dji_backward
clear d2i d2j

%% -------------------------------------------------------------- ALGORITHM
% (1) INITIALIZATION: u^0 = 0, z^0 = 0, solve for c=1,2,...
u   = reshape(input,[],C);
chi = reshape(mask,[],C);

param.lambda = param.lambda*chi;

rhsL = (param.lambda/param.gamma).*reshape(u,[],C);
rhsM = ones(M*N,C); 

% FOR EACH COLOR CHANNEL
for c = 1:C
    
    % ITERATION
    for iter = 1:maxiter
      
        % SOLVE EULER-LAGRANGE EQUATION FOR \chi
        % i.e M(u^{c-1},\chi^c) = 1.
        % M is a linear operator acting on \chi and reads
        % M(u,.) = 1+\frac{2\epsilon\alpha}{\beta} |\nabla u|^2
        %           - 4\epsilon^2\Delta.
        % Solved via inversion of the linear operators.
        MM       = matrixM(param,M*N,matrices,u(:,c));
        chinew   = MM\rhsM(:,c);
        diff_chi = norm(chinew-chi(:,c))/norm(chinew);
        chi(:,c) = chinew;
        clear chinew
        
        % SOLVE EULER-LAGRANGE EQUATION FOR u: 
        % i.e. L(\chi^c,u^c)= \alpha \chi_{\Omega\setminus D} \cdot ustart.
        % L is a linear operator acting on u and reads
        % L(\chi,.) = -\div(\chi_\epsilon^2\nabla)
        %             + \alpha\chi_{\Omega\setminus D}.
        % Solved via inversion of the linear operators.
        LL     = matrixL(param,M*N,matrices,chi(:,c));
        unew   = LL\rhsL(:,c);
        diff_u = norm(unew-u(:,c))/norm(unew);
        u(:,c) = unew;
        clear unew
        
        % TEST EXIT CONDITION
        if diff_u<tol
            break
        end    
        
    end    
end

u   = reshape(u,M,N,C);
chi = reshape(chi,M,N,C);

%% WRITE IMAGE OUTPUTS
imwrite(u,'./results/mumford_shah_output.png')
%imwrite(chi,'./results/mumford_shah_levels_output.png')

return

%% ---------------------------------------------------- AUXILIARY FUNCTIONS
function M = matrixM(param,N,matrices,u)
% Definition of (\nabla u)^2:
nablau2 = (matrices.Dic*u).^2 + (matrices.Djc*u).^2;

M = speye(N)...
    + 2 * param.epsilon * param.gamma/param.alpha * spdiags(nablau2,0,N,N)...
    - 4*param.epsilon^2*matrices.LAP;
return

function L = matrixL(param,N,matrices,chi)
% Definition of the nonlinear diffusion weighted by \chi^2:

z  = chi.^2 + param.epsilon^2; % coefficient of nonlinear diffusion

zx = matrices.Dic*z;
zy = matrices.Djc*z;

Z  = spdiags(z, 0,N,N);
Zx = spdiags(zx,0,N,N);
Zy = spdiags(zy,0,N,N);

NonlinearDelta = Z*matrices.LAP + Zx*matrices.Dic + Zy*matrices.Djc;

L =  -NonlinearDelta + spdiags(param.lambda/param.gamma,0,N,N);

return
