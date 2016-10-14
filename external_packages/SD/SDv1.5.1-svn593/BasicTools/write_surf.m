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
%
% write_surf(vertex_coords, faces, fname)
%
% Assume vertex_coords is vnum x 3
% writes the vertex coordinates and face lists of a triangular mesh to a
% surface file
%

function write_surf(vertex_coords, faces, fname)



faces = int32(faces);
TRIANGLE_FILE_MAGIC_NUMBER =  16777214 ;

fid = fopen(fname, 'wb', 'b') ;
if(fid < 0)
    error('could not open surface file %s to write', fname) ;
end

fwrite3(fid, TRIANGLE_FILE_MAGIC_NUMBER);

fprintf(fid, '\n');
fprintf(fid, '\n');
 

vnum = int32(size(vertex_coords, 1));
fnum = int32(size(faces, 1));

fwrite(fid, vnum, 'int32');
fwrite(fid, fnum, 'int32');

fwrite(fid, vertex_coords', 'float32');
fwrite(fid, faces', 'int32');
fclose(fid);
