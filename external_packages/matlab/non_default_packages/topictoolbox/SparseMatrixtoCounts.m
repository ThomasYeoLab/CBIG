function [ WS , DS ] = SparseMatrixtoCounts( WD )
%% Function SparseMatrixtoCounts

%%
% [ WS , DS ] = SparseMatrixtoCounts( WD )

%%
% 

[ ii , jj , ss ] = find( WD );

ntokens = full( sum( sum( WD )));
WS = zeros( 1,ntokens );
DS = zeros( 1,ntokens );

startindex = 0;
for i=1:length( ii )
   nt = ss( i );
   WS( startindex+1:startindex+nt ) = ii( i );
   DS( startindex+1:startindex+nt ) = jj( i );
   startindex = startindex + nt; 
end


