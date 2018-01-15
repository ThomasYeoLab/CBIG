function b = subsref(ipd,index)
% INPLACE/SUBSREF Support Matlab subscript references.
%
% This function supports subscript references just like Matlab matrices.
%
% Example:
%    ipa = inplace(cumsum(ones(5,1)));
%    ipa([2 4])
%    ipa(3:end)
    
if (length(index) > 1)
    error('Indexing operation not supported');
end;
if (strcmp(index.type,'()'))
    b = subsref(ipd.get_a(),index);
elseif (strcmp(index.type,'.'))
    try
        a = 1;
        a.a;
    catch
        rethrow(lasterr);
    end;
elseif (strcmp(index.type,'{}'))
    try
        a = 1;
        a{1};
    catch
        rethrow(lasterr);
    end;
else
    error('Unsupported');
end;
