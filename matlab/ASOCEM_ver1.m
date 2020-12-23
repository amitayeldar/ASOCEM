
function ASOCEM_ver1(micrograph_addr,output_dir,particle_size,downscale_size,area_size,smoothing_term)
% parameters
files = dir([micrograph_addr,'/','*.mrc']);
numOfMicro = size(files,1);
% numOfMicro = 1;
maxIterAsocem=600;
parfor expNum = 1:numOfMicro
    [~, microName] = fileparts(files(expNum).name);
    mgBig = ReadMRC([files(expNum).folder,'/',files(expNum).name]);
    mgBig = double(mgBig);
    mgBig = rot90(mgBig);
    I0 = mgBig;
    % run ASOCEM
    [phi] = ASOCEM(I0,downscale_size,area_size,smoothing_term,maxIterAsocem);
    if phi==ones(size(phi)) % all is contamination dont pick
        continue
    end
    % get rid of blubs the size of the particle
    phi_seg = imbinarize(zeros(size(phi)));
    scalingSz = downscale_size/max(size(I0));
    min_bulb_size = floor(2*scalingSz*particle_size/2);
    % se_dil = strel('disk',ceil(max(1,scalingSz*particle_size/8)));
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
    end 
    % resizing phi_seg to original image
    phi_seg = imresize(phi_seg,size(mgBig));
    f=figure('visible', 'off');
    subplot(1,2,1); imshow(cryo_downsample(I0,200),[]);
    subplot(1,2,2); imshow(imresize(phi_seg,[200,200]),[]);
    saveas(f,[output_dir,microName,'_imageSz_',num2str(downscale_size),'area_',num2str(area_size),'.jpg'])
%     save([output_dir,microName,'_binary_seg.mat'],'phi_seg');
end
end