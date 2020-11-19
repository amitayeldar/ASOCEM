%% script for ASOCEM
% ASOCEM = automated segmentation of contamination elctron microscopy.

% clear all
addpath('./ASOCEM_functions');
%% real data sets
% loading the micrograph
img_adrr = './Data/h200/GridSquare_9221822_Data_FoilHole_9224247_Data_9224878_9224879_20181109_1515-1447_aligned_mic_DW.mrc';
I0 = ReadMRC(img_adrr);
I0 = double(I0);

% parameters
band_pass = 1;
area_mat_sz = 9; % after down scaling to 200*200
smoothing_term = 1-band_pass; % a number from 0 to 1. 0 not to smooth and 1 to max smoothing
maxIter=200;
particle_size = 200;

% run ASOCEM
[phi] = ASOCEM(I0,area_mat_sz,smoothing_term,maxIter);

% get rid of blubs the size of the particle
phi_seg = zeros(size(phi));
scalingSz = 200/max(size(I0));
min_bulb_size = floor(2*scalingSz*particle_size/2);
se_dil = strel('disk',1);
se_erod = strel('disk',min_bulb_size);
CC = bwconncomp(phi>0,8);
for i =1:size(CC.PixelIdxList,2)
       if size(CC.PixelIdxList{i},1)> (scalingSz*particle_size)^2 
                tmp = zeros(size(phi));
                tmp(CC.PixelIdxList{i})= 1;
                tmp_dil = imdilate(tmp,se_dil);
                tmp_erode = imerode(tmp_dil,se_erod);
                if nnz(tmp_erode(:)) > 0
                    phi_seg(CC.PixelIdxList{i})=1;
                end
       end
        t = ceil(area_mat_sz/2+1);
        % we do not know about the boundary of the image hence is contamination
        phi_seg(1:t,:) = 1;
        phi_seg(end-t:end,:) = 1;
        phi_seg(:,1:t) = 1;
        phi_seg(:,end-t:end) = 1;
        % resizing phi_seg to original image
        phi_seg = imresize(phi_seg,size(I0));
end



%% resizing phi to original image
% phi = cryo_downsample(phi,size(I0));
% phi_seg = cryo_downsample(phi_seg,size(I0));



%% figures
figure;suptitle(['areaMatSz = ',num2str(area_mat_sz),' bandPass ',num2str(band_pass)]);
subplot(1,2,1); imshow(cryo_downsample(I0,size(phi)),[])
subplot(1,2,2); imshow(phi_seg,[])
saveas(gcf,[img_adrr(1:end-4),'ams',num2str(area_mat_sz),'.png'],'png')

%% correlated variable with radial covarience
% 
% % parameters
% area_mat_sz = 11; % odd number only
% if mod(area_mat_sz,2)==0 %% cov_mat_sz have to be odd so cov matrix size will be odd
%     area_mat_sz = area_mat_sz-1;
% end
% initMethod = 0; % 0 for starting with circ, 1 for starting with gauss log test
% smoothing = 0;
% 
% maxImgSz = 201;
% maxIter=200;
% tol= 0.001;
% 
% cov_mat_sz = area_mat_sz^2;
% max_d = floor(sqrt(2)*area_mat_sz)+1;
% 
% % creating the image
% % background
% max_d = ceil(sqrt(2)*9);
% sz1 = 199;
% mu1 = 1;
% rad_deacey = 10; % the exp coeff in noise_exp2d
% u1 = noise_exp2d_amitay(sz1,1,rad_deacey);
% [R1_true,x,cnt]=cryo_epsdR(u1,1:sz1^2,max_d); % computing real radial_cov;
% u1 =u1 - mean(u1(:)) + mu1;
% 
% % contamination
% sz0 = 199;
% mu0 = 0.5;
% rad_deacey = 0.01; % the exp coeff in noise_exp2d
% u0 = noise_exp2d_amitay(sz0,1,rad_deacey);
% [R0_true,x,cnt]=cryo_epsdR(u0,1:sz0^2,max_d);
% u0 =u0 - mean(u0(:)) + mu0;
% % 
% % % prelocating contaminations
% phi_true = zeros(sz1);
% [X,Y] = meshgrid(1:sz1);
% cx = 60; cy = 100; %center of circ
% rowIdxB = X-cx;
% colIdxB = Y-cy;
% Rsquare = (X-cx).^2+(Y-cy).^2;
% phi_true(Rsquare<=(30^2)) = 1;
% phi_true(70:100,20:90) = 1;
% phi_true(10:70,170:185) = 1;
% cx = 190; cy = 190; %center of circ
% rowIdxB = X-cx;
% colIdxB = Y-cy;
% Rsquare = (X-cx).^2+(Y-cy).^2;
% phi_true(Rsquare<=(30^2)) = 1;
% 
% u1(phi_true==1) = 0;
% u0(phi_true==0) = 0;
% I0 = u1+u0;
% 
% 
% % run chan vese
% 
% [I0,phi0,phi,mu0_est,R0_est,mu1_est,R1_est] = segmentation(I0,cov_mat_sz,maxImgSz,initMethod,smoothing,maxIter,tol);
% 
% 
% % err
% err_R0 = norm(R0_true-R0_est)/norm(R0_true);
% err_mu0 = abs(mu0-mu0_est)/abs(mu0);
% err_R1 = norm(R1_true-R1_est)/norm(R1_true);
% err_mu1 = abs(mu1-mu1_est)/abs(mu1);
% phi_tmp = phi>0;
% phi_true_tmp = phi_seg;
% phi_diff = abs(phi_seg-phi_true);
% err_phi = sum(phi_diff(:))/sum(phi_true(:));
% 
% 
% % figures
% figure;suptitle(['areaMatSz = ',num2str(area_mat_sz),' meanSep = ',num2str(abs(mu0-mu1))]);
% subplot(2,2,1); imshow(I0,[])
% subplot(2,2,2); imshow(phi_true>0)
% subplot(2,2,3); imshow(phi0>0)
% subplot(2,2,4); imshow(phi>0)


%% independet variable
% % sz1 = 201;
% % s1 = 10;
% % mu1 = 1*ones(sz1^2,1);
% % cov1 = diag(sqrt(s1)*ones(sz1^2,1));
% % u1 = cov1 * randn(sz1^2,1) + mu1;
% % u1 = reshape(u1,[sz1,sz1]);
% % sz0 = 50;
% % s0 = 1;
% % mu0 = 1*ones(sz0^2,1);
% % cov0 = diag(sqrt(s0)*ones(sz0^2,1));
% % u0 = cov0 * randn(sz0^2,1) + mu0;
% % u0 = reshape(u0,[sz0,sz0]);
% % 
% % 
% % 
% % I0 = u1;
% % I0(10:10+sz0-1,10:10+sz0-1) = u0;
% % I0(100:100+sz0-1,100:100+sz0-1) = u0;
% % 
% % % region = zeros(201);
% % % region(10:10+sz0-1,10:10+sz0-1) = 1;
% % % region(100:100+sz0-1,100:100+sz0-1) = 1;
% % % sample_idx = find(region>0);
% % % 
% % %  [R,x,cnt]=cryo_epsdR(I0-20,sample_idx,max_d);
% % 
% % [phi,I0,phi0,mu0_est,s0_est,mu1_est,s1_est] = segmentation(I0,-1,max_d,1,0,0.001,0);
% % %% err
% % err_mu0 = abs(mean(mu0)-mu0_est)/abs(mean(mu0))
% % err_s0 = abs(s0-mean(diag(s0_est)))/abs(s0)
% % err_mu1 = abs(mean(mu1)-mu1_est)/abs(mean(mu1))
% % err_s1 = abs(s1-mean(diag(s1_est)))/abs(s1)
% % 
% % % figures
% % figure; subplot(1,3,1); imshow(I0,[])
% % subplot(1,3,2); imshow(phi0>0)
% % subplot(1,3,3); imshow(phi>0)
% 
% 
% %% real exp
% clear all
% 
% max_d  = 6;
% img_adrr  = './Data/006';
% img_adrr = [img_adrr,'.mrc'];
% [phi,I0,phi0,mu0_est,R0_est,mu1_est,R1_est] = segmentation(img_adrr,-1,max_d,1,0,0.001,0);
% 
% % figures
% figure;suptitle(['maxd = ',num2str(max_d)]);
% subplot(1,3,1); imshow(I0,[])
% subplot(1,3,2); imshow(phi0>0)
% subplot(1,3,3); imshow(phi>0)
% 

