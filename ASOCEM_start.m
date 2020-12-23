function ASOCEM_start

% ASOCEM_start
% Gather all information required to start the contamination removal out of the micrographs.
% This is the first command to be called in any processing workflw.
% 
% Amitay Eldar, December 2020.

if ~isdeployed % Only run in a MATLAB session
    [basedir,~,~]=fileparts(mfilename('fullpath'));
    addpath(fullfile(basedir,'matlab')); % set up MATLAB path
end

micrograph_addr='';
while isempty(micrograph_addr)
    micrograph_addr =fmtinput('Enter full path of micrographs MRC file: ','','%s');
%     if exist(micrograph_addr,'file')~=7
    if isempty(dir([micrograph_addr,'/*.mrc']))
        fprintf('MRC file does not exist.\n');
        micrograph_addr='';
    end
end

output_dir =fmtinput('Enter full path of output directory: ','','%s');
if ~strcmp(output_dir(end),'/')
    output_dir = [output_dir,'/'];
end
    
if ~exist(output_dir,'dir') % Do we need to create directory?
    message='Output directory does not exist. Create?';
    do_create=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
    if do_create==1
        mkdir(output_dir);
    end
end

particle_size='';
while isempty(particle_size)
    particle_size_str =fmtinput('Enter the particle size in pixels: ','','%s');
    particle_size = str2double(particle_size_str);
    if mod(particle_size,1)~=0
        fprintf('particle size should be a natural number.\n');
        particle_size='';
    end
    if particle_size<0
        fprintf('particle size should be a natural number.\n');
        particle_size='';
    end
end

downscale_size='';
while isempty(downscale_size)
    downscale_size_str =fmtinput('Enter the image size after downscaling in pixels: ','','%s');
    downscale_size = str2double(downscale_size_str);
    if mod(downscale_size,1)~=0
        fprintf('downscale size should be a natural number.\n');
        downscale_size='';
    end
    if downscale_size<0
        fprintf('downscale size should be a natural number.\n');
        downscale_size='';
    end
end

area_size='';
while isempty(area_size)
    area_size_str =fmtinput('Enter the area size after downscaling in pixels: ','','%s');
    area_size = str2double(area_size_str);
    if mod(area_size,1)~=0
        fprintf('Area size should be a natural number.\n');
        area_size='';
    end
    if area_size<0
        fprintf('Area size should be a natural number.\n');
        area_size='';
    end
end


smoothing_term='';
while isempty(smoothing_term)
    smoothing_term_str =fmtinput('Enter the smoothing_term 1 for max smoothing and 0 not to smooth: ','','%s');
    smoothing_term = str2double(smoothing_term_str);
    if or(smoothing_term<0,smoothing_term>1)
        fprintf('smoothing_term should be a number between 0 to 1.\n');
        smoothing_term='';
    end
end

% message='Do you want to use the GPU?';
% do_create=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
% if do_create==1
%     gpu_use = 1;
% else
%     gpu_use = 0;
% end


ASOCEM_ver1(micrograph_addr,output_dir,particle_size,downscale_size,area_size,smoothing_term)

