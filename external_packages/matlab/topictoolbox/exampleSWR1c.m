%% Example 1 of running special word retrieval model (SWR)
%

%%
% Choose the dataset
dataset = 1; % 1 = psych review abstracts 2 = NIPS papers

if (dataset == 1)
    % Set the hyperparameters 
    GAMMA = [ 0.01 0.01 0.01 ];
    
    BURNIN = 20;  % the number of iterations
    NS     = 100;  % number of samples
    LAG    = 1; % lag between samples
    
    SEED   = 1;  % The random seed
    OUTPUT = 0; % What output to show (0=no output; 1=iterations; 2=all output)

    % load the psych review data in bag of words format  
    load 'psychreviewstream'; % DS WS WO
    WOS = WO;
    DSS = DS;
    WSS = WS;
    load 'psychreviewbagofwordsfromstream'; % WO DS WS ORIGWSPOS;  
    load 'swrsamples_1\\swrdata'; % PQ0_ALLW PQ1_ALLW PQ2_ALLW
    D = size( PQ0_ALLW , 2 );
    for i=1:D
        titles{ i } = sprintf( 'Document %d\n' , i );
    end
elseif (dataset == 2)
    % Set the hyperparameters 
    GAMMA = [ 0.01 0.01 0.01 ]; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    BURNIN = 20;  % the number of iterations
    NS     = 100;  % number of samples
    LAG    = 1; % lag between samples
    
    SEED   = 1;  % The random seed
    OUTPUT = 0; % What output to show (0=no output; 1=iterations; 2=all output)

    % load the psych review data in bag of words format 
    load 'nips_stream'; % DS WS WO
    WOS = WO;
    DSS = DS;
    WSS = WS;
    load 'nipsbagofwordsfromstream'; % WO DS WS ORIGWSPOS;  
    load 'swrsamples_2\\swrdata'; % PQ0_ALLW PQ1_ALLW PQ2_ALLW
    %load titles_nips;
    D = size( PQ0_ALLW , 2 );
    for i=1:D
        titles{ i } = sprintf( 'Document %d\n' , i );
    end
end

% create document word sparse count matrix
WD = sparse( WS , DS , ones( size( DS )));
W = size( WD , 1 );
D = size( WD , 2 );

ok = 1;

while (ok==1)
    inputquery = input( 'Query: ' , 's' );
    inputquery = deblank( lower( inputquery ));
    
    if length( inputquery ) == 0
        ok = 0;
    else
        Query = tokenize( inputquery );

        % translate query into indices
        Query = lower( Query );
        NQ = length( Query );
        Q = zeros( 1,NQ );
        for i=1:NQ
            wh = find( strcmp( Query{ i } , WO ));
            if isempty( wh )
                fprintf( 'Unknown word: %s\n' , Query{ i } );
                wh = 0;
            end
            Q( i ) = wh;
        end

        if all( Q > 0 )

            nmatch = full( WD( Q , : ));

            whallmatch = sum( nmatch > 0 , 1 ) == NQ;
            fprintf( '\nRetrieval results: there are %d docs with all keywords\n\n' , sum( whallmatch ));
            
            % create the NQ x D word probabilities according to topic model
            PQ0 = PQ0_ALLW( Q , : );
            PQ1 = PQ1_ALLW( Q , : );
            PQ2 = PQ2_ALLW( Q );

            % Run the mex routine
            [ X , XD , PQ ] = GibbsSamplerSWR( PQ0 , PQ1 , PQ2 , GAMMA , BURNIN , NS , LAG , SEED , OUTPUT );

%             Weights = prod( PQ , 1 );
%             X2 = zeros( 3 , NQ );
%             if NQ > 1
%                 for d=1:D
%                     X2 = X2 + Weights( d ) * squeeze( XD( d , : , : ));
%                 end
%             else
%                 for d=1:D
%                     X2 = X2 + Weights( d ) * squeeze( XD( d , :  ))';
%                 end
%             end
            
%             fprintf( 'Overall route probabilities [topic,doc,corpus]:\n' );
%             for i=1:NQ
%                 fprintf( '%20s ' , Query{ i } );
%                 weights = X2( : , i );
%                 weights = weights / sum( weights );
%                 for r=1:3
%                     fprintf( '%4.3f  ' , weights( r ));
%                 end
%                 fprintf( '\n' );
%             end
%             fprintf( '\n' );

            Score = sum( -log( PQ ) , 1 );
            [ Score , index ] = sort( Score );
            
            d = 0;
            moreresults = 1;
            while (d < D) & (moreresults==1)
                d = d + 1;
            
                doc = index( d );
                fprintf( 'Rank:%2d Score=%4.4f  ' , d , Score( d ) );
                fprintf( 'Whole query routes=[%4.3f %4.3f %4.3f]\n\n' , X( doc , 2 ) , X( doc , 1 ) , X( doc , 3 ));
                for i=1:NQ
                    fprintf( '\t%14s Nmatches=%2d prob=%4.3f routes=[%4.3f %4.3f %4.3f]\n' , sprintf( '"%s"' , Query{ i }) , nmatch( i , doc ) , PQ( i , doc ) , XD( doc , 2 , i ) , XD( doc , 1 , i ) , XD( doc , 3 , i ));
                end
                fprintf( '\n' );
                fprintf( '%s\n' , titles{ doc } );
                
                % show a fragment of the abstract
                ii = 0;
                whind = find( DSS == doc );
                wsind = WSS( whind );
                nwdoc = length( wsind );
                maxlength = 250;
                ll = 0;
                cap = 1;
                while (ii < min( nwdoc , maxlength))
                    ii = ii + 1;
                    if wsind( ii ) == 0
                        str = '.';
                        cap = 1;
                    else
                        str = lower( WOS{ wsind( ii )});
                        if any( strcmp( str , Query ))
                            str = [ '*' upper( str ) '*' ];
                        end
                        if (cap==1)
                            str(1) = upper( str(1));
                        end
                        cap = 0;
                    end
                    ll = ll + length( str ) + 1;
                    
                    if wsind( ii ) == 0
                       fprintf( '%s' , str );
                    else
                       fprintf( ' %s' , str ); 
                    end
                    
                    
                    if ll > 80
                        ll = 0;
                        fprintf( '\n' );
                    end
                end
                
                fprintf( '\n\n' );
                
                inputquery = input( 'More search results? [enter=YES; any letter=NO]: ' , 's' );
                if length( inputquery ) > 0
                    moreresults = 0;
                end
                
                fprintf( '\n' );
            end
        else
           fprintf( '\n' ); 
        end
    end
end





