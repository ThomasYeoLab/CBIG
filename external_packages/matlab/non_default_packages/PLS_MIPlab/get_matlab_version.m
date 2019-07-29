function vernum = get_matlab_version

   ab = version;
   dot = findstr(ab, '.');
   a = ab(1:dot(1)-1);
   b = ab(dot(1)+1:dot(2)-1);
   vernum = str2num(a)*1000+str2num(b);

   return;

