% Contact ythomas@csail.mit.edu or msabuncu@csail.mit.edu for bugs or questions 
%
%=========================================================================
%
%  Copyright (c) 2008 Thomas Yeo and Mert Sabuncu
%  All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%
%    * Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
%
%    * Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
%    * Neither the names of the copyright holders nor the names of future
%      contributors may be used to endorse or promote products derived from this
%      software without specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
%ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.    
%
%=========================================================================
function zi = MARS_interp2(varargin)
%INTERP2 2-D interpolation (table lookup).
%   ZI = INTERP2(X,Y,Z,XI,YI) interpolates to find ZI, the values of the
%   underlying 2-D function Z at the points in matrices XI and YI.
%   Matrices X and Y specify the points at which the data Z is given.
%
%   XI can be a row vector, in which case it specifies a matrix with
%   constant columns. Similarly, YI can be a column vector and it 
%   specifies a matrix with constant rows. 
%
%   ZI = INTERP2(Z,XI,YI) assumes X=1:N and Y=1:M where [M,N]=SIZE(Z).
%   ZI = INTERP2(Z,NTIMES) expands Z by interleaving interpolates between
%   every element, working recursively for NTIMES.  INTERP2(Z) is the
%   same as INTERP2(Z,1).
%
%   ZI = INTERP2(...,METHOD) specifies alternate methods.  The default
%   is linear interpolation.  Available methods are:
%
%     'nearest' - nearest neighbor interpolation
%     'linear'  - bilinear interpolation
%     'cubic'   - bicubic interpolation
%     'spline'  - spline interpolation
%
%   For faster interpolation when X and Y are equally spaced and monotonic,
%   use the syntax ZI = INTERP2(...,*METHOD).
%
%   ZI = INTERP2(...,METHOD,EXTRAPVAL) specificies a method and a scalar 
%   value for ZI outside of the domain created by X and Y.  Thus, ZI will
%   equal EXTRAPVAL for any value of YI or XI which is not spanned by Y 
%   or X respectively. A method must be specified for EXTRAPVAL to be used,
%   the default method is 'linear'.
%
%   All the interpolation methods require that X and Y be monotonic and
%   plaid (as if they were created using MESHGRID).  If you provide two
%   monotonic vectors, interp2 changes them to a plaid internally. 
%   X and Y can be non-uniformly spaced.
%
%   For example, to generate a coarse approximation of PEAKS and
%   interpolate over a finer mesh:
%       [x,y,z] = peaks(10); [xi,yi] = meshgrid(-3:.1:3,-3:.1:3);
%       zi = interp2(x,y,z,xi,yi); mesh(xi,yi,zi)
%
%   Class support for inputs X, Y, Z, XI, YI:  
%      float: double, single
%
%   See also INTERP1, INTERP3, INTERPN, MESHGRID, GRIDDATA.

%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 5.33.4.11 $
error(nargchk(1,7,nargin)); % allowing for an ExtrapVal

bypass = false;
uniform = true;
if (nargin > 1)
    if nargin == 7 && ~isnumeric(varargin{end})
        error('MATLAB:interp2:extrapvalNotNumeric',...
              'EXTRAPVAL must be numeric.');
    end
    if ischar(varargin{end})
        narg = nargin-1;
        method = [varargin{end} '    ']; % Protect against short string.
        ExtrapVal = nan; % setting default ExtrapVal as NAN
        index = 1; %subtract off the elements not in method
    elseif ischar(varargin{end-1}) && isnumeric(varargin{end})
        narg = nargin-2;
        method = [ varargin{end-1} '    '];
        ExtrapVal = varargin{end}; % user specified ExtrapVal
        index = 2; % subtract off the elements not in method and ExtrapVal
    else
        narg = nargin;
        method = 'linear';
        ExtrapVal = nan; % protecting default
        index = 0;
    end
    if strncmpi(method,'*',1) % Direct call bypass.
        if (narg ==5 || narg ==3)
            xitemp = varargin{end-index - 1};
            yitemp = varargin{end-index};
            if isrowvector(xitemp) && iscolvector(yitemp)
                varargin{end-index - 1} = repmat(xitemp, [size(yitemp,1), 1]);
                varargin{end-index} = repmat(yitemp, [1, size(xitemp,2)]);
            elseif iscolvector(xitemp) && isrowvector(yitemp)
                varargin{end-index - 1} = repmat(xitemp, [1, size(yitemp, 2)]);
                varargin{end-index} = repmat(yitemp, [size(xitemp,1), 1]);
            end
        end
        if strcmpi(method(2),'l') || strcmpi(method(2:4),'bil') 
            % bilinear interpolation.
            zi = linear(ExtrapVal, varargin{1:end-index});
            return
        elseif strcmpi(method(2),'c') || strcmpi(method(2:4),'bic') 
            % bicubic interpolation
            zi = cubic(ExtrapVal, varargin{1:end-index});
            return
        elseif strcmpi(method(2),'n') 
            % Nearest neighbor interpolation
            zi = nearest(ExtrapVal, varargin{1:end-index});
            return
        elseif strcmpi(method(2),'s') 
            % spline interpolation
            method = 'spline'; bypass = true;
        else
            error('MATLAB:interp2:InvalidMethod',...
                  [deblank(method),' is an invalid method.']);
        end
    elseif strncmpi(method,'s',1), % Spline interpolation
        method = 'spline'; bypass = true;
    end
else
    narg = nargin;
    method = 'linear';
    ExtrapVal = nan; % default ExtrapVal is NaN
end

if narg==1, % interp2(z), % Expand Z
    [nrows,ncols] = size(varargin{1});
    xi = 1:.5:ncols; yi = (1:.5:nrows)';
    x = 1:ncols; y = 1:nrows;
    [msg,x,y,z,xi,yi] = xyzchk(x,y,varargin{1},xi,yi);

elseif narg==2. % interp2(z,n), Expand Z n times
    [nrows,ncols] = size(varargin{1});
    ntimes = floor(varargin{2}(1));
    xi = 1:1/(2^ntimes):ncols; yi = (1:1/(2^ntimes):nrows)';
    x = 1:ncols; y = 1:nrows;
    [msg,x,y,z,xi,yi] = xyzchk(x,y,varargin{1},xi,yi);

elseif narg==3, % interp2(z,xi,yi)
    [nrows,ncols] = size(varargin{1});
    x = 1:ncols; y = 1:nrows;
    [msg,x,y,z,xi,yi] = xyzchk(x,y,varargin{1:3});

elseif narg==4,
    error('MATLAB:interp2:nargin','Wrong number of input arguments.');

elseif narg==5, % linear(x,y,z,xi,yi)
    [msg,x,y,z,xi,yi] = xyzchk(varargin{1:5});

end

error(msg);

%
% Check for plaid data.
%
xx = x(1,:); yy = y(:,1);
if (size(x,2)>1 && ~isequal(repmat(xx,size(x,1),1),x)) || ...
   (size(y,1)>1 && ~isequal(repmat(yy,1,size(y,2)),y)),
    error('MATLAB:interp2:meshgrid',...
        ['X and Y must be matrices produced by MESHGRID. Use' ...
        ' GRIDDATA instead \nof INTERP2 for scattered data.']);
end

%
% Check for non-equally spaced data.  If so, map (x,y) and
% (xi,yi) to matrix (row,col) coordinate system.
%
if ~bypass,
    xx = xx.'; % Make sure it's a column.
    dx = diff(xx); dy = diff(yy);
    xdiff = max(abs(diff(dx))); if isempty(xdiff), xdiff = 0; end
    ydiff = max(abs(diff(dy))); if isempty(ydiff), ydiff = 0; end
    if (xdiff > eps*max(abs(xx))) || (ydiff > eps*max(abs(yy))),
        if any(dx < 0), % Flip orientation of data so x is increasing.
            x = fliplr(x); y = fliplr(y); z = fliplr(z);
            xx = flipud(xx); dx = -flipud(dx);
        end
        if any(dy < 0), % Flip orientation of data so y is increasing.
            x = flipud(x); y = flipud(y); z = flipud(z);
            yy = flipud(yy); dy = -flipud(dy);
        end

        if any(dx<=0) || any(dy<=0),
            error('MATLAB:interp2:XorYNotMonotonic',...
                  'X and Y must be monotonic vectors or matrices produced by MESHGRID.');
        end

        % Bypass mapping code for cubic
        if ~strncmp(method(1),'c',1)
            % Determine the nearest location of xi in x
            [xxi,j] = sort(xi(:));
            [ignore,i] = sort([xx;xxi]);
            ui(i) = 1:length(i);
            ui = (ui(length(xx)+1:end)-(1:length(xxi)))';
            ui(j) = ui;

            % Map values in xi to index offset (ui) via linear interpolation
            ui(ui<1) = 1;
            ui(ui>length(xx)-1) = length(xx)-1;
            ui = ui + (xi(:)-xx(ui))./(xx(ui+1)-xx(ui));

            % Determine the nearest location of yi in y
            [yyi,j] = sort(yi(:));
            [ignore,i] = sort([yy;yyi(:)]);
            vi(i) = 1:length(i);
            vi = (vi(length(yy)+1:end)-(1:length(yyi)))';
            vi(j) = vi;

            % Map values in yi to index offset (vi) via linear interpolation
            vi(vi<1) = 1;
            vi(vi>length(yy)-1) = length(yy)-1;
            vi = vi + (yi(:)-yy(vi))./(yy(vi+1)-yy(vi));

            [x,y] = meshgrid(ones(class(x)):size(x,2),ones(class(y)):size(y,1));
            xi(:) = ui; yi(:) = vi;
        else
            uniform = false;
        end
    end
end

% Now do the interpolation based on method.
if strncmpi(method,'l',1) || strncmpi(method,'bil',3) % bilinear interpolation.
    zi = linear(ExtrapVal,x,y,z,xi,yi);

elseif strncmpi(method,'c',1) || strncmpi(method,'bic',3) % bicubic interpolation
    if uniform
        zi = cubic(ExtrapVal,x,y,z,xi,yi);
    else
        zi = spline2(x,y,z,xi,yi,ExtrapVal);
        d = xi < min(x(:)) | xi > max(x(:)) | yi < min(y(:)) | yi > max(y(:));        
        zi(d) = ExtrapVal;
    end

elseif strncmpi(method,'n',1) % Nearest neighbor interpolation
    zi = nearest(ExtrapVal,x,y,z,xi,yi);

elseif strncmpi(method,'s',1) % Spline interpolation
    zi = spline2(x,y,z,xi,yi,ExtrapVal);

else
    error('MATLAB:interp2:InvalidMethod',...
        [deblank(method),' is an invalid method.']);

end

%------------------------------------------------------
function irv = isrowvector(rv)
% check if input is row vector
irv = isvector(rv) && size(rv,1)==1;

%------------------------------------------------------
function icv = iscolvector(rv)
% check if input is column vector
icv = isvector(rv) && size(rv,2)==1;

%------------------------------------------------------
function F = linear(ExtrapVal,arg1,arg2,arg3,arg4,arg5)
%LINEAR 2-D bilinear data interpolation.
%   ZI = LINEAR(EXTRAPVAL,X,Y,Z,XI,YI) uses bilinear interpolation to
%   find ZI, the values of the underlying 2-D function in Z at the points
%   in matrices XI and YI.  Matrices X and Y specify the points at which
%   the data Z is given.  X and Y can also be vectors specifying the
%   abscissae for the matrix Z as for MESHGRID. In both cases, X
%   and Y must be equally spaced and monotonic.
%
%   Values of EXTRAPVAL are returned in ZI for values of XI and YI that are
%   outside of the range of X and Y.
%
%   If XI and YI are vectors, LINEAR returns vector ZI containing
%   the interpolated values at the corresponding points (XI,YI).
%
%   ZI = LINEAR(EXTRAPVAL,Z,XI,YI) assumes X = 1:N and Y = 1:M, where
%   [M,N] = SIZE(Z).
%
%   ZI = LINEAR(EXTRAPVAL,Z,NTIMES) returns the matrix Z expanded by
%   interleaving bilinear interpolates between every element, working
%   recursively for NTIMES. LINEAR(EXTRAPVAL,Z) is the same as
%   LINEAR(EXTRAPVAL,Z,1).
%
%   See also INTERP2, CUBIC.

if nargin==2 % linear(extrapval,z), Expand Z
    [nrows,ncols] = size(arg1);
    s = 1:.5:ncols; lengths = length(s);
    t = (1:.5:nrows)'; lengtht = length(t);
    s = repmat(s,lengtht,1);
    t = repmat(t,1,lengths);
    
elseif nargin==3 % linear(extrapval,z,n), Expand Z n times
    [nrows,ncols] = size(arg1);
    ntimes = floor(arg2);
    s = 1:1/(2^ntimes):ncols; lengths = length(s);
    t = (1:1/(2^ntimes):nrows)'; lengtht = length(t);
    s = repmat(s,lengtht,1);
    t = repmat(t,1,lengths);

elseif nargin==4 % linear(extrapval,z,s,t), No X or Y specified.
    [nrows,ncols] = size(arg1);
    s = arg2; t = arg3;

elseif nargin==5
    error('MATLAB:interp2:linear:nargin','Wrong number of input arguments.');

elseif nargin==6 % linear(extrapval,x,y,z,s,t), X and Y specified.
    [nrows,ncols] = size(arg3);
    mx = numel(arg1); my = numel(arg2);
    if (mx ~= ncols || my ~= nrows) && ~isequal(size(arg1),size(arg2),size(arg3))
        error('MATLAB:interp2:linear:XYZLengthMismatch',...
            'The lengths of the X and Y vectors must match Z.');
    end
    if nrows < 2 || ncols < 2
        error('MATLAB:interp2:linear:sizeZ','Z must be at least 2-by-2.');
    end
    s = 1 + (arg4-arg1(1))/(arg1(end)-arg1(1))*(ncols-1);
    t = 1 + (arg5-arg2(1))/(arg2(end)-arg2(1))*(nrows-1);

end

if nrows < 2 || ncols < 2
    error('MATLAB:interp2:linear:sizeZsq','Z must be at least 2-by-2.');
end
if ~isequal(size(s),size(t))
    error('MATLAB:interp2:linear:XIandYISizeMismatch',...
        'XI and YI must be the same size.');
end

% Check for out of range values of s and set to 1
sout = find((s<1)|(s>ncols));
if length(sout)>0, s(sout) = 1; end

% Check for out of range values of t and set to 1
tout = find((t<1)|(t>nrows));
if length(tout)>0, t(tout) = 1; end

% Matrix element indexing
ndx = floor(t)+floor(s-1)*nrows;

% Compute intepolation parameters, check for boundary value.
if isempty(s), d = s; else d = find(s==ncols); end
s(:) = (s - floor(s));
if length(d)>0, s(d) = s(d)+1; ndx(d) = ndx(d)-nrows; end

% Compute intepolation parameters, check for boundary value.
if isempty(t), d = t; else d = find(t==nrows); end
t(:) = (t - floor(t));
if length(d)>0, t(d) = t(d)+1; ndx(d) = ndx(d)-1; end

% Now interpolate.
onemt = 1-t;
if nargin==6,
    F =  ( arg3(ndx).*(onemt) + arg3(ndx+1).*t ).*(1-s) + ...
         ( arg3(ndx+nrows).*(onemt) + arg3(ndx+(nrows+1)).*t ).*s;
else
    F =  ( arg1(ndx).*(onemt) + arg1(ndx+1).*t ).*(1-s) + ...
         ( arg1(ndx+nrows).*(onemt) + arg1(ndx+(nrows+1)).*t ).*s;
end

% Now set out of range values to NaN.
if length(sout)>0, F(sout) = ExtrapVal; end
if length(tout)>0, F(tout) = ExtrapVal; end

%------------------------------------------------------
function F = cubic(ExtrapVal,arg1,arg2,arg3,arg4,arg5)
%CUBIC 2-D bicubic data interpolation.
%   CUBIC(...) is the same as LINEAR(....) except that it uses
%   bicubic interpolation.
%
%   This function needs about 7-8 times SIZE(XI) memory to be available.
%
%   See also LINEAR.

%   Based on "Cubic Convolution Interpolation for Digital Image
%   Processing", Robert G. Keys, IEEE Trans. on Acoustics, Speech, and
%   Signal Processing, Vol. 29, No. 6, Dec. 1981, pp. 1153-1160.

if nargin==2, % cubic(extrapval,z), Expand Z
    [nrows,ncols] = size(arg1);
    s = 1:.5:ncols; lengths = length(s);
    t = (1:.5:nrows)'; lengtht = length(t);
    s = repmat(s,lengtht,1);
    t = repmat(t,1,lengths);
    
elseif nargin==3, % cubic(extrapval,z,n), Expand Z n times
    [nrows,ncols] = size(arg1);
    ntimes = floor(arg2);
    s = 1:1/(2^ntimes):ncols; lengths = length(s);
    t = (1:1/(2^ntimes):nrows)'; lengtht = length(t);
    s = repmat(s,lengtht,1);
    t = repmat(t,1,lengths);

elseif nargin==4, % cubic(extrapval,z,s,t), No X or Y specified.
    [nrows,ncols] = size(arg1);
    s = arg2; t = arg3;

elseif nargin==5,
    error('MATLAB:interp2:cubic:nargin','Wrong number of input arguments.');

elseif nargin==6, % cubic(extrapval,x,y,z,s,t), X and Y specified.
    [nrows,ncols] = size(arg3);
    mx = numel(arg1); my = numel(arg2);
    if (mx ~= ncols || my ~= nrows) && ~isequal(size(arg1),size(arg2),size(arg3))
        error('MATLAB:interp2:cubic:XYZLengthMismatch',...
            'The lengths of the X and Y vectors must match Z.');
    end
    if nrows < 3 || ncols < 3
        error('MATLAB:interp2:cubic:sizeZ','Z must be at least 3-by-3.');
    end
    s = 1 + (arg4-arg1(1))/(arg1(end)-arg1(1))*(ncols-1);
    t = 1 + (arg5-arg2(1))/(arg2(end)-arg2(1))*(nrows-1);

end

if nrows < 3 || ncols < 3
    error('MATLAB:interp2:cubic:sizeZsq','Z must be at least 3-by-3.');
end
if ~isequal(size(s),size(t)),
    error('MATLAB:interp2:cubic:XIandYISizeMismatch','XI and YI must be the same size.');
end

% Check for out of range values of s and set to 1
sout = find((s<1)|(s>ncols));
if length(sout)>0, s(sout) = 1; end

% Check for out of range values of t and set to 1
tout = find((t<1)|(t>nrows));
if length(tout)>0, t(tout) = 1; end

% Matrix element indexing
ndx = floor(t)+floor(s-1)*(nrows+2);

% Compute intepolation parameters, check for boundary value.
if isempty(s), d = s; else d = find(s==ncols); end
s(:) = (s - floor(s));
if length(d)>0, s(d) = s(d)+1; ndx(d) = ndx(d)-nrows-2; end

% Compute intepolation parameters, check for boundary value.
if isempty(t), d = t; else d = find(t==nrows); end
t(:) = (t - floor(t));
if length(d)>0, t(d) = t(d)+1; ndx(d) = ndx(d)-1; end

if nargin==6,
    % Expand z so interpolation is valid at the boundaries.
    zz = zeros(size(arg3)+2);
    zz(1,2:ncols+1) = 3*arg3(1,:)-3*arg3(2,:)+arg3(3,:);
    zz(2:nrows+1,2:ncols+1) = arg3;
    zz(nrows+2,2:ncols+1) = 3*arg3(nrows,:)-3*arg3(nrows-1,:)+arg3(nrows-2,:);
    zz(:,1) = 3*zz(:,2)-3*zz(:,3)+zz(:,4);
    zz(:,ncols+2) = 3*zz(:,ncols+1)-3*zz(:,ncols)+zz(:,ncols-1);
    nrows = nrows+2; %also ncols = ncols+2;
else
    % Expand z so interpolation is valid at the boundaries.
    zz = zeros(size(arg1)+2);
    zz(1,2:ncols+1) = 3*arg1(1,:)-3*arg1(2,:)+arg1(3,:);
    zz(2:nrows+1,2:ncols+1) = arg1;
    zz(nrows+2,2:ncols+1) = 3*arg1(nrows,:)-3*arg1(nrows-1,:)+arg1(nrows-2,:);
    zz(:,1) = 3*zz(:,2)-3*zz(:,3)+zz(:,4);
    zz(:,ncols+2) = 3*zz(:,ncols+1)-3*zz(:,ncols)+zz(:,ncols-1);
    nrows = nrows+2; %also ncols = ncols+2;
end

% Now interpolate using computationally efficient algorithm.
t0 = ((2-t).*t-1).*t;
t1 = (3*t-5).*t.*t+2;
t2 = ((4-3*t).*t+1).*t;
t(:) = (t-1).*t.*t;
F     = ( zz(ndx).*t0 + zz(ndx+1).*t1 + zz(ndx+2).*t2 + zz(ndx+3).*t ) ...
    .* (((2-s).*s-1).*s);
ndx(:) = ndx + nrows;
F(:)  = F + ( zz(ndx).*t0 + zz(ndx+1).*t1 + zz(ndx+2).*t2 + zz(ndx+3).*t ) ...
    .* ((3*s-5).*s.*s+2);
ndx(:) = ndx + nrows;
F(:)  = F + ( zz(ndx).*t0 + zz(ndx+1).*t1 + zz(ndx+2).*t2 + zz(ndx+3).*t ) ...
    .* (((4-3*s).*s+1).*s);
ndx(:) = ndx + nrows;
F(:)  = F + ( zz(ndx).*t0 + zz(ndx+1).*t1 + zz(ndx+2).*t2 + zz(ndx+3).*t ) ...
    .* ((s-1).*s.*s);
F(:) = F/4;

% Now set out of range values to NaN.
if length(sout)>0, F(sout) = ExtrapVal; end
if length(tout)>0, F(tout) = ExtrapVal; end

%------------------------------------------------------
function F = nearest(ExtrapVal,arg1,arg2,arg3,arg4,arg5)
%NEAREST 2-D Nearest neighbor interpolation.
%   ZI = NEAREST(EXTRAPVAL,X,Y,Z,XI,YI) uses nearest neighbor interpolation
%   to find ZI, the values of the underlying 2-D function in Z at the points
%   in matrices XI and YI.  Matrices X and Y specify the points at which
%   the data Z is given.  X and Y can also be vectors specifying the
%   abscissae for the matrix Z as for MESHGRID. In both cases, X
%   and Y must be equally spaced and monotonic.
%
%   Values of EXTRAPVAL are returned in ZI for values of XI and YI that are
%   outside of the range of X and Y.
%
%   If XI and YI are vectors, NEAREST returns vector ZI containing
%   the interpolated values at the corresponding points (XI,YI).
%
%   ZI = NEAREST(EXTRAPVAL,Z,XI,YI) assumes X = 1:N and Y = 1:M, where
%   [M,N] = SIZE(Z).
%
%   F = NEAREST(EXTRAPVAL,Z,NTIMES) returns the matrix Z expanded by
%   interleaving interpolates between every element.  NEAREST(EXTRAPVAL,Z)
%   is the same as NEAREST(EXTRAPVAL,Z,1).
%
%   See also INTERP2, LINEAR, CUBIC.

if nargin==2, % nearest(z), Expand Z
    [nrows,ncols] = size(arg1);
    u = 1:.5:ncols; lengthu = length(u);
    v = (1:.5:nrows)'; lengthv = length(v);
    u = repmat(u,lengthv,1);
    v = repmat(v,1,lengthu);

elseif nargin==3, % nearest(z,n), Expand Z n times
    [nrows,ncols] = size(arg1);
    ntimes = floor(arg2);
    u = 1:1/(2^ntimes):ncols; lengthu = length(u);
    v = (1:1/(2^ntimes):nrows)'; lengthv = length(v);
    u = repmat(u,lengthv,1);
    v = repmat(v,1,lengthu);

elseif nargin==4, % nearest(z,u,v)
    [nrows,ncols] = size(arg1);
    u = arg2; v = arg3;

elseif nargin==5,
    error('MATLAB:interp2:nearest:nargin','Wrong number of input arguments.');

elseif nargin==6, % nearest(x,y,z,u,v), X and Y specified.
    [nrows,ncols] = size(arg3);
    mx = numel(arg1); my = numel(arg2);
    if (mx ~= ncols || my ~= nrows) && ...
            ~isequal(size(arg1),size(arg2),size(arg3))
        error('MATLAB:interp2:nearest:XYZLengthMismatch',...
            'The lengths of the X and Y vectors must match Z.');
    end
    if nrows > 1 && ncols > 1
        u = 1 + (arg4-arg1(1))/(arg1(mx)-arg1(1))*(ncols-1);
        v = 1 + (arg5-arg2(1))/(arg2(my)-arg2(1))*(nrows-1);
    else
        u = 1 + (arg4-arg1(1));
        v = 1 + (arg5-arg2(1));
    end
end

if ~isequal(size(u),size(v))
    error('MATLAB:interp2:nearest:XIandYISizeMismatch',...
        'XI and YI must be the same size.');
end

% Check for out of range values of u and set to 1
uout = (u<.5)|(u>=ncols+.5);
anyuout = any(uout(:));
if anyuout, u(uout) = 1; end

% Check for out of range values of v and set to 1
vout = (v<.5)|(v>=nrows+.5);
anyvout = any(vout(:));
if anyvout, v(vout) = 1; end

% Interpolation parameters
u = round(u); v = round(v);

% Now interpolate
ndx = v+(u-1)*nrows;
if nargin==6,
    F = arg3(ndx);
else
    F = arg1(ndx);
end

% Now set out of range values to NaN.
if anyuout, F(uout) = ExtrapVal; end
if anyvout, F(vout) = ExtrapVal; end

%----------------------------------------------------------
function F = spline2(varargin)
%2-D spline interpolation

% Determine abscissa vectors
varargin{1} = varargin{1}(1,:);
varargin{2} = varargin{2}(:,1).';

%
% Check for plaid data.
%
xi = varargin{4}; yi = varargin{5};
xxi = xi(1,:); yyi = yi(:,1);

if ~isequal(repmat(xxi,size(xi,1),1),xi) || ...
        ~isequal(repmat(yyi,1,size(yi,2)),yi)
    F = splncore(varargin(2:-1:1),varargin{3},varargin(5:-1:4));
else
    F = splncore(varargin(2:-1:1),varargin{3},{yyi(:).' xxi},'gridded');
end
