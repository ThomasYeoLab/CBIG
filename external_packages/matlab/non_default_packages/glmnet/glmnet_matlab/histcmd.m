function cmd = histcmd()

fid = fopen([prefdir,'/history.m'],'rt');
while ~feof(fid)
   last = fgetl(fid);
end
fclose(fid);
cmd = last;

end