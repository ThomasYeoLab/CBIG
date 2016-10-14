%% Example 3 for utilizing the LDA model
% This example shows how to order topics according to co-usage in documents
% This will an ordering over topics where topics that co-appear often with
% other topics in documents to have similar indices.
%
% This examples uses the OrderTopics function. In this procedure, for each
% topic, the probability distribution over documents is calculated. For
% each topic pair (i,j), the (symmetrized) KL distance between topic
% distributions i and j is calculated. The KL distances are then subjected
% multidimensional scaling (MDS). The first dimension of the MDS solution
% is then used to get an ordering for topics.
% 
% load a document-topic count matrix saved for the nips dataset
load 'ldasingle_nips';
load 'words_nips';

%%
% extract the topics in a cell array of strings
[S]=WriteTopics( WP,BETA,WO,5,0.6 );

%%
% order these topics using the OrderTopics function
[ Order ] = OrderTopics( DP,ALPHA );

%%
% and show the resulting ordering
S( Order )