% script to create contaminated micrograph with radial auto corr for 
% both background and contamination
clear all
if ~isdeployed % Only run in a MATLAB session
    [basedir,~,~]=fileparts(mfilename('fullpath'));
    addpath(fullfile(basedir,'functions')); % set up MATLAB path
end

%% correlated variable with general radial covarience
% parameters
imgSz = 400;
mu0mean = 1;
var0 =1 ;
a0 = 0.1;
mu1mean = 4;
var1 = 1;
a1 = 1;

num_of_sqr_cont = 1;
min_max_sz_sqr = [0.1,0.2]*imgSz;
num_of_circ_cont = 1;
min_max_sz_circ = [0.1,0.2]*imgSz;

path_to_save_image = basedir;
image_name = ['/CVWGRC_patchSz_',num2str(imgSz),'.mrc'];

% creating distance matrix
cov_mat_sz = imgSz^2;
X = zeros(cov_mat_sz,2);
cnt = 1;
for i=1:sqrt(cov_mat_sz)
    for j=1:sqrt(cov_mat_sz)
        X(cnt,1) = i;
        X(cnt,2) = j;
        cnt = cnt + 1;
    end
end

Z = squareform(pdist(single(X)));

% creating u0,u1 statistics
E0 = exp(-a0*Z);
E1 = exp(-a1*Z);

mu0 = normrnd(mu0mean,0.1,cov_mat_sz,1);
cov0 = E0-diag(diag(E0)-var0);
mu1 = normrnd(mu1mean,0.1,cov_mat_sz,1);
cov1 = E1-diag(diag(E1)-var1);

% creating u0,u1
u0 = reshape(mvnrnd(mu0,cov0,1),imgSz,imgSz);
u1 = reshape(mvnrnd(mu1,cov1,1),imgSz,imgSz);

% % prelocating contaminations
phi_true = zeros(imgSz);
[X,Y] = meshgrid(1:imgSz);
sqr_location = randi([1,imgSz],num_of_sqr_cont,2);
sqr_sz = randi(min_max_sz_sqr,num_of_sqr_cont,1);
circ_location = randi([1,imgSz],num_of_circ_cont,2);
circ_sz = randi(min_max_sz_circ,num_of_circ_cont,1);
for i=1:num_of_sqr_cont
    Xtmp = abs(X-X(sqr_location(i,1),sqr_location(i,2)));
    Xtmp = Xtmp<sqr_sz(i,1);
    Ytmp = abs(Y-Y(sqr_location(i,1),sqr_location(i,2)));
    Ytmp = Ytmp<sqr_sz(i,1);
    phi_true(Xtmp.*Ytmp==1) = 1;
end

for i=1:num_of_circ_cont
    Xtmp = X-X(circ_location(i,1),circ_location(i,2));
    Ytmp = Y-Y(circ_location(i,1),circ_location(i,2));
    phi_true(Xtmp.^2+Ytmp.^2<circ_sz(i,1)^2) = 1;
end
u0tmp = u0; u1tmp=u1;
u0tmp(phi_true<=0) = 0;
u1tmp(phi_true>0) = 0;
I0 = u0tmp + u1tmp;
figure; imshow(I0,[])
WriteMRC(I0,1,[path_to_save_image,image_name]);


























%% correlated variable with independables patches

% parameters
imgSz = 1000;
patchSz = 3;
num_of_patch = floor(imgSz/patchSz);
mu0mean = 1;
var0mean = 1;
mu1mean = 4;
var1mean = 1;

num_of_sqr_cont = 1;
min_max_sz_sqr = [0.1,0.2]*imgSz;
num_of_circ_cont = 1;
min_max_sz_circ = [0.1,0.2]*imgSz;

path_to_save_image = basedir;
image_name = ['/CVWIP_patchSz_',num2str(patchSz),'.mrc'];

% making u0 and u1
% u0 parameters
mu0 = normrnd(mu0mean,var0mean,patchSz^2,1)';
cov0 = (1/50)*randi([-50,50],patchSz^2,patchSz^2);
cov0 = 0.5*(cov0+cov0');
cov0 = cov0+3*var0mean*eye(patchSz^2);
% u1 parameters
mu1 = normrnd(mu1mean,var1mean,patchSz^2,1)';
cov1 = (1/50)*randi([-50,50],patchSz^2,patchSz^2);
cov1 = 0.5*(cov1+cov1');
cov1 = cov1+3*var1mean*eye(patchSz^2);
u0=[];u1=[];
for i=1:num_of_patch
    tmp0=[];tmp1=[];
    for j=1:num_of_patch
        tmp0 = [tmp0,reshape(mvnrnd(mu0,cov0,1),patchSz,patchSz)];
        tmp1 = [tmp1,reshape(mvnrnd(mu1,cov1,1),patchSz,patchSz)];
    end
    if i==1
        u0=tmp0;
        u1=tmp1;
    else
        u0=[u0;tmp0];
        u1=[u1;tmp1];    
    end
end
% % prelocating contaminations
img_sz=size(u0,1);
phi_true = zeros(img_sz);
[X,Y] = meshgrid(1:img_sz);
sqr_location = randi([1,img_sz],num_of_sqr_cont,2);
sqr_sz = randi(min_max_sz_sqr,num_of_sqr_cont,1);
circ_location = randi([1,img_sz],num_of_circ_cont,2);
circ_sz = randi(min_max_sz_circ,num_of_circ_cont,1);
for i=1:num_of_sqr_cont
    Xtmp = abs(X-X(sqr_location(i,1),sqr_location(i,2)));
    Xtmp = Xtmp<sqr_sz(i,1);
    Ytmp = abs(Y-Y(sqr_location(i,1),sqr_location(i,2)));
    Ytmp = Ytmp<sqr_sz(i,1);
    phi_true(Xtmp.*Ytmp==1) = 1;
end

for i=1:num_of_circ_cont
    Xtmp = X-X(circ_location(i,1),circ_location(i,2));
    Ytmp = Y-Y(circ_location(i,1),circ_location(i,2));
    phi_true(Xtmp.^2+Ytmp.^2<circ_sz(i,1)^2) = 1;
end

u0(phi_true<=0) = 0;
u1(phi_true>0) = 0;
I0 = u0 + u1;
figure; imshow(I0,[])
WriteMRC(I0,1,[path_to_save_image,image_name]);



%% correlated variable with radial covarience of yoel
path_to_save_image = basedir;
image_name = '/syn_mic.mrc';
% parameters
img_sz = 800;
mu0 = 0.5;
var0 = 1;
rad_deacey_coeff_0 = 1;
mu1 = 4;
var1 = 1;
rad_deacey_coeff_1 = 10;


num_of_sqr_cont = 1;
min_max_sz_sqr = [0.1,0.2]*img_sz;
num_of_circ_cont = 1;
min_max_sz_circ = [0.1,0.2]*img_sz;

% creating the image
% contamination
u0 = noise_exp2d_amitay(img_sz,1,rad_deacey_coeff_0,sqrt(var0));
u0 =u0 - mean(u0(:)) + mu0;
% background
u1 = noise_exp2d_amitay(img_sz,1,rad_deacey_coeff_1,sqrt(var1));
u1 =u1 - mean(u1(:)) + mu1;
% % prelocating contaminations
phi_true = zeros(img_sz);
[X,Y] = meshgrid(1:img_sz);
sqr_location = randi([1,img_sz],num_of_sqr_cont,2);
sqr_sz = randi(min_max_sz_sqr,num_of_sqr_cont,1);
circ_location = randi([1,img_sz],num_of_circ_cont,2);
circ_sz = randi(min_max_sz_circ,num_of_circ_cont,1);
for i=1:num_of_sqr_cont
    Xtmp = abs(X-X(sqr_location(i,1),sqr_location(i,2)));
    Xtmp = Xtmp<sqr_sz(i,1);
    Ytmp = abs(Y-Y(sqr_location(i,1),sqr_location(i,2)));
    Ytmp = Ytmp<sqr_sz(i,1);
    phi_true(Xtmp.*Ytmp==1) = 1;
end

for i=1:num_of_circ_cont
    Xtmp = X-X(circ_location(i,1),circ_location(i,2));
    Ytmp = Y-Y(circ_location(i,1),circ_location(i,2));
    phi_true(Xtmp.^2+Ytmp.^2<circ_sz(i,1)^2) = 1;
end

u0(phi_true<=0) = 0;
u1(phi_true>0) = 0;
I0 = u0 + u1;
figure; imshow(I0,[])
WriteMRC(I0,1,[path_to_save_image,image_name]);
