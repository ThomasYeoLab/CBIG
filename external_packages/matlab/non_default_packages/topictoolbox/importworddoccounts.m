function [ WS , DS ] = importworddoccounts( filename , nheaderlines , str );
%% Imports text file with word-document counts into matlab format
%
%%
%
% [ WS , DS ] = importworddoccounts( filename ) reads in word counts per documents and
% produces the word and document indices |WS| and |DS| respectively. |WS(k)| and |DS(k)| contains the kth word and
% document index respectively. Note that each repetition of a word in a
% document will be represented by a separate token in the |WS| and |DS|
% vectors.
%
%%
%
% [ WS , DS ] = importworddoccounts( filename , nheaderlines ) will skip the first
% nheaderlines lines of the file
%
% [ WS , DS ] = importworddoccounts( filename , nheaderlines , 'reverse' ) will
% reverse the word and documents columns
%
%%
% NOTES
% The text file should be organized into
% three colums where each row contains the document index, the word index,
% and the word count. For example:
%      1 2 10
%      1 3 4 
%      2 2 6
% should be read as "word 2 occurs 10 times in doc 1, word 3 occurs 4 times
% in doc 1, and word 2 occurs 6 times in doc 2".

if nargin==1
    nheaderlines = 0;
end

doreverse = 0;
if (nargin==3) & strcmp( str , 'reverse' )
    doreverse = 1;
end
     
fprintf( 'Importing word document counts from: %s\n' , filename );

if (doreverse==0)
    [ dd,ww,cc ] = textread( filename , '%d %d %d' , 'headerlines' , nheaderlines );
    WD = sparse( ww,dd,cc );
else
    fprintf( 'Reversing word-document columns\n' );
    [ ww,dd,cc ] = textread( filename , '%d %d %d' , 'headerlines' , nheaderlines );
    WD = sparse( ww,dd,cc );
end

fprintf( 'Number of Documents D = %d\n' , size( WD , 2 ));
fprintf( 'Number of Words     W = %d\n' , size( WD , 1 ));
fprintf( 'Number of nonzero entries in matrix NNZ=%d\n' , nnz( WD ));

[ WS , DS ] = SparseMatrixtoCounts( WD );

