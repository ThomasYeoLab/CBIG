function ip = inplace(a)
% INPLACE Convert to a type that supports inplace modification
% ipa = inplace(v) creates an inplace double object from any matrix v.
%
% Example:
%    ipa = inplace(ones(5));

if nargin == 0
    error('inplace:invalidParameter', ['ipdouble must be created with an initial .']);
elseif isa(a,'ipdouble')
    % make a copy
    ipd = inplace(a.a);
    
    return;
    
end

    function out = ip_get_a()
        out = a;
    end

    function ip_assign(in)
        a = in;
    end

    function ip_subsasgn(S,B)
        subsasgn(a,S,B);
    end

    ip.get_a = @ip_get_a;
    ip.assign = @ip_assign;
    ip.subsasgn = @ip_subsasgn;
    
    ip = class(ip,'inplace');

% ***** end inplace ***** 
end



