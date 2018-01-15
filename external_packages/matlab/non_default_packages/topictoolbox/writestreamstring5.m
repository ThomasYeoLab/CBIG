function writestreamstring5( WS , WO , ncharsperline , maxlines , fid , x , z , thr , NC , T );
  
L = [];
N = length( WS );
NL = 0;
for i=1:N
    if (WS( i ) == 0)
        L      = [ L '. ' ];
        spanst = 'normal';
        word   = '. ';
    else
        word   = WO{ WS( i )};
        probs  = x( : , i );
        topic  = z( i );
        
        [ temp , maxi ] = max( probs );
        
        L      = [ L ' ' word ];
                
        if (topic==0)
            spanst = 'filler'; 
        elseif (maxi==1) & (topic <= T)
            spanst = 'normal';
            word = [ word sprintf( '<sup>%d</sup>' , topic ) ]; 
        elseif (maxi==2) | (topic == T+1)
            spanst = 'normal';
            whcol = round( probs( 2 ) * NC );
            spanst = sprintf( 't%d' , whcol );
        else
            %word = [ word sprintf( '<sup>*</sup>' )]; 
            spanst = 'backg';
        end
    end
    
    if length( L ) > ncharsperline
        NL = NL + 1;
        
        if (NL == maxlines+1)
            fprintf( fid , '...' );
            break
        else
            L = word;
            fprintf( fid , '\r\n' );
        end
    end

    fprintf( fid , ' <span class="%s">%s</span>' , spanst , word );
end

fprintf( fid , '\r\n' ); 

