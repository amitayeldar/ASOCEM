function [phi,mu0_est,cov0_est,mu1_est,cov1_est] = ASOCEM(I0,downscale_size,area_mat_sz,smoothing_term,maxIter)
% This function segment contamination in cryo_EM micrographs

% INPUT
% I0               micrograph image .mrc format
% area_mat_sz      area for covarience estimation   
% smoothing_term   a number between 0 to 1. 0 not to smooth and 1 for
%                  maximum smoothing
% maxIter          maximum number of iteration allowed

% OUTPUT
% phi          chan vese level set function. phi>0 is area zero and phi<=0 is area 1
% mu0_est      area zero mean estimation
% R0_est       area zero radial covarience estimation
% mu1_est      area one mean estimation
% R1_est       area one radial covarience estimation

%% intialize segmentation parameters
dt = 10^(0); % time step
mu = 10^(-2); %curve length regularization
nu = 0; %curve area regularization
eta = 10^(-8); %prevent zero curvature (division by 0)
Eps = 1; %epsilon (for delta function)
tol = 10^(-3);
% maxImgSz = 200; % downsampling parameter
if mod(area_mat_sz,2)==0 %% cov_mat_sz have to be odd so cov matrix size will be odd
    area_mat_sz = area_mat_sz-1;
end
cov_mat_sz = area_mat_sz^2;

%% preprocess image
% size rescale
szScaling = downscale_size/max(size(I0));
I = cryo_downsample(I0,floor(szScaling*size(I0)));
% image size should be odd
if mod(size(I,1),2)==0
    I = I(1:end-1,:);
end
if mod(size(I,2),2)==0
    I = I(:,1:end-1);
end

% smoothing
if smoothing_term>0 && smoothing_term<1
    if 1-smoothing_term <= 0.3
         bandPass1d = fir1(10, [0.001 0.3]);
         bandPass2d  = ftrans2(bandPass1d); %% radial bandpass
         I = imfilter(I, bandPass2d,'replicate');
    else
         bandPass1d = fir1(10, [0.001 1-smoothing_term]);
         bandPass2d  = ftrans2(bandPass1d); %% radial bandpass
         I = imfilter(I, bandPass2d,'replicate');
    end
end

%% intialize phi as liphsiczh circ
[X,Y] = meshgrid(-floor(size(I,2)/2):1:floor(size(I,2)/2),-floor(size(I,1)/2):1:floor(size(I,1)/2));
% phi_0 = sin((pi/30)*X).*sin((pi/30)*Y);
phi_0 = min(size(I)/3)^2 - (X.^2 + Y.^2);
phi_0 = phi_0./max(abs(phi_0(:)));


%% chan vese time process
[phi,mu0_est,cov0_est,mu1_est,cov1_est] = chan_vese_process(I,phi_0,cov_mat_sz,dt,mu,nu,eta,Eps,maxIter,tol);
if phi == ones(size(phi)) % we dont want to use this micrograph
    return
end

% the smaller area will considered to be the contamination

size_0 = sum(phi>0,'all');
size_1 = sum(phi<0,'all');
if size_0 >= size_1 % changing areas 
    phi = -1*phi;
end
  
end

