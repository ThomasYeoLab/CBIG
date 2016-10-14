dataset = 2;

if (dataset==1)
    stopfile  = '..\kdd\NYTimes\stopwordlistNYTimes.txt';
    load 'nips_stream'; % WS DS WO;
    [ WW , WO , WS , DS , SI , ORIGWSPOS ]=stream_to_collocation_data(DS,WS,WO,stopfile);
    save 'nipsbagofwordsfromstream' WO DS WS ORIGWSPOS;
elseif (dataset==2)
    stopfile  = '..\kdd\NYTimes\stopwordlistNYTimes.txt';
    load 'psychreviewstream'; % WS DS WO;
    [ WW , WO , WS , DS , SI , ORIGWSPOS ]=stream_to_collocation_data(DS,WS,WO,stopfile);
    save 'psychreviewbagofwordsfromstream' WO DS WS ORIGWSPOS;
end
