function [varargout] = rdir(rootdir,varargin)
% Lists the files in a directory and its sub directories. 
% 
% function [D] = rdir(ROOT,TEST)
%
% Recursive directory listing.
%
% ROOT is the directory starting point and includes the 
% wildcard specification.
% The function returns a structure D similar to the one 
% returned by the built-in dir command. 
% There is one exception, the name field will include 
% the relative path as well as the name to the file that 
% was found.
% Pathnames and wildcards may be used. Wild cards can exist
% in the pathname too. A special case is the double * that
% will match multiple directory levels, e.g. c:\**\*.m. 
% Otherwise a single * will only match one directory level.
% e.g. C:\Program Files\Windows *\
%
% TEST is an optional test that can be performed on the 
% files. Two variables are supported, datenum & bytes.
% Tests are strings similar to what one would use in a 
% if statement. e.g. 'bytes>1024 & datenum>now-7'
%
% If not output variables are specified then the output is 
% sent to the screen.
%
% See also DIR
%
% examples:
%   D = rdir('*.m');
%     for ii=1:length(D), disp(D(ii).name); end;
%
%   % to find all files in the current directory and sub directories
%   D = rdir('**\*')
%
%   % If no output is specified then the files are sent to 
%   % the screen.
%   rdir('c:\program files\windows *\*.exe');
%   rdir('c:\program files\windows *\**\*.dll');
%
%   % Using the test function to find files modified today
%   rdir('c:\win*\*','datenum>floor(now)');
%   % Using the test function to find files of a certain size
%   rdir('c:\program files\win*\*.exe','bytes>1024 & bytes<1048576');
%

% use the current directory if nothing is specified
if ~exist('rootdir','var'),
  rootdir = '*';
end;

% split the file path around the wild card specifiers
prepath = '';       % the path before the wild card
wildpath = '';      % the path wild card
postpath = rootdir; % the path after the wild card
I = find(rootdir==filesep,1,'last');
if ~isempty(I),
  prepath = rootdir(1:I);
  postpath = rootdir(I+1:end);
  I = find(prepath=='*',1,'first');
  if ~isempty(I),
    postpath = [prepath(I:end) postpath];
    prepath = prepath(1:I-1);
    I = find(prepath==filesep,1,'last');
    if ~isempty(I),
      wildpath = prepath(I+1:end);
      prepath = prepath(1:I);
    end;
    I = find(postpath==filesep,1,'first');
    if ~isempty(I),
      wildpath = [wildpath postpath(1:I-1)];
      postpath = postpath(I:end);
    end;
  end;
end;
% disp([' "' prepath '" ~ "' wildpath '" ~ "' postpath '" ']);


if isempty(wildpath),
  % if no directory wildcards then just get file list
  D = dir([prepath postpath]);
  D([D.isdir]==1) = [];
  for ii = 1:length(D),
    if (~D(ii).isdir),
      D(ii).name = [prepath D(ii).name];
    end;
  end;

  % disp(sprintf('Scanning "%s"   %g files found',[prepath postpath],length(D)));

elseif strcmp(wildpath,'**'), % a double wild directory means recurs down into sub directories

  % first look for files in the current directory (remove extra filesep)
  D = rdir([prepath postpath(2:end)]);

  % then look for sub directories
  Dt = dir(''); 
  tmp = dir([prepath '*']);
  % process each directory
  for ii = 1:length(tmp),
    if (tmp(ii).isdir && ~strcmpi(tmp(ii).name,'.') && ~strcmpi(tmp(ii).name,'..') ),
      Dt = [Dt; rdir([prepath tmp(ii).name filesep wildpath postpath])];
    end;
  end;
  D = [D; Dt];

else
  % Process directory wild card looking for sub directories that match
  tmp = dir([prepath wildpath]);
  D = dir(''); 
  % process each directory found
  for ii = 1:length(tmp),
    if (tmp(ii).isdir && ~strcmpi(tmp(ii).name,'.') && ~strcmpi(tmp(ii).name,'..') ),
      D = [D; rdir([prepath tmp(ii).name postpath])];
    end;
  end;
end;


% Apply filter
if (nargin>=2 && ~isempty(varargin{1})),
  date = [D.date];
  datenum = [D.datenum];
  bytes = [D.bytes];

  try
    eval(sprintf('D((%s)==0) = [];',varargin{1})); 
  catch
    warning('Error: Invalid TEST "%s"',varargin{1});
  end;
end;

% display listing if no output variables are specified
if nargout==0,
  pp = {'' 'k' 'M' 'G' 'T'};
  for ii=1:length(D), 
    sz = D(ii).bytes;
    if sz<=0,
      disp(sprintf(' %31s %-64s','',D(ii).name)); 
    else
      ss = min(4,floor(log2(sz)/10));
      disp(sprintf('%4.0f %1sb   %20s   %-64s ',sz/1024^ss,pp{ss+1},D(ii).date,D(ii).name)); 
    end;
  end;
else
  % send list out
  varargout{1} = D;
end;

