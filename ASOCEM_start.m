function ASOCEM_start

% ASOCEM_start
% Gather all information required to start the contamination removal out of the micrographs.
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

contamination_criterion='';
while isempty(contamination_criterion)
    contamination_criterion_str =fmtinput('Enter the contamination_criterion, 0 for lower size, 1 to lower mean: ','','%s');
    contamination_criterion = str2double(contamination_criterion_str);
    if and(contamination_criterion~=0,contamination_criterion~=1)
        fprintf('contamination_criterion should be  0 or 1.\n');
        contamination_criterion='';
    end
end

fast_flag='';
while isempty(fast_flag)
    fast_flag_str =fmtinput('Enter flag_fast 0,1,2 where 2 is the fastest: ','','%s');
    fast_flag = str2double(fast_flag_str);
    if and(and(fast_flag~=0,fast_flag~=1),fast_flag~=2)
        fprintf('fast_flag should be  0,1,2 .\n');
        fast_flag='';
    end
end



ASOCEM_ver1(micrograph_addr,output_dir,particle_size,downscale_size,area_size,contamination_criterion,fast_flag)

