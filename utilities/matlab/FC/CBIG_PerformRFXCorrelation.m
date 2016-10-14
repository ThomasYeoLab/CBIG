function CBIG_PerformRFXCorrelation(pfile, input_file_txt, exit_flag)

% CBIG_PerformRFXCorrelation(pfile, input_file_txt, exit_flag)
%
% This function is used to compute p-val of correlation greater than 0
% across data on input_file_txt
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(nargin < 3)
   exit_flag = 0; 
else
   if(ischar(exit_flag))
      exit_flag = str2num(exit_flag); 
   end
end

fid = fopen(input_file_txt, 'r');
i = 0;
while(1);
    tmp = fscanf(fid, '%s\n', 1);
    if(isempty(tmp))
        break
    else
        i = i + 1;
        input_files{i} = tmp;
    end
end
fclose(fid);

for i = 1:length(input_files)
    
    x = MRIread(input_files{i});
    correlation = reshape(x.vol, size(x.vol, 1)*size(x.vol,2)*size(x.vol,3), size(x.vol, 4));
    
    tmp = CBIG_StableAtanh(correlation); % fisher-z transform
    if(i == 1)
        zmat = tmp;
        zsqmat = tmp.^2;
    else
        zmat = zmat + tmp; %update sum x
        zsqmat = zsqmat + tmp.^2; % update sum x^2
    end
end
zmat = zmat/length(input_files); % compute mean
zsqmat = zsqmat/length(input_files);
stdz = sqrt(zsqmat - zmat.^2)*length(input_files)/(length(input_files) - 1); %compute std
tmat = sqrt(length(input_files))*zmat./stdz; % compute t stats
pmat = tcdf(-abs(tmat), length(input_files) - 1); % compute p val

pout = x;
pout.vol = reshape(pmat, size(x.vol));
MRIwrite(pout, pfile);

if(exit_flag)
   exit; 
end
