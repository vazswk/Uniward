%% mian.m
clc;clear all;
currdir = pwd;
indir=[currdir '\images_cover\'];
outdir=[currdir '\images_stego\'];
%mkdir(outdir);
payload = single(0.2);

flist = dir([indir '\*.jpg']);
flen = length(flist);
fprintf('%s%d\n', 'the num of the files: ',flen);

% if exist(outdir,'dir'); rmdir(outdir,'s'); end    
% if ~exist(outdir,'dir'); mkdir(outdir); end
   for i = 1:flen
       fprintf('%d%s\n',i, ['   processing image: ' flist(i).name]);
       cover_path = [indir flist(i).name];
       stego_path = [outdir flist(i).name];

       config.STC_h = uint32(10);
       config.seed  = int32(123);
       distortion = J_UNIWARD(cover_path, stego_path, payload,config);
      
        
   end

 
 