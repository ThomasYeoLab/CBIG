%% DATA REQUIREMENTS
%
% To get started on applying topic models to your own corpus of text, there
% are various Matlab variables that need to be constructed. The exact format
% depends on which type of topic model you want to run. Below is a list of
% specification for the input that needs to be provided to each kind of
% topic model and the kind of output each model will produce
%
%% LDA MODEL
%
%%
% The input is a bag of word representation containing the number of times each
% words occurs in a document. The outputs are the topic assignments to each
% word token as well as the counts of the number of times each word is
% assigned to a topic and the number of times a topic is assigned to a
% document
%
%%
% *INPUT*
% 
% * *|WS|*  a 1 x |N| vector where |WS(k)| contains the vocabulary index of
% the kth word token, and |N| is the number of word tokens. The word
% indices are not zero based, i.e., min( |WS| )=1 and max( |WS| ) = |W| =
% number of distinct words in vocabulary
%
% * *|DS|*  a 1 x |N| vector where |DS(k)| contains the document index of
% the kth word token. The document indices are not zero based, i.e., min(
% |DS| )=1 and max( |DS| ) = |D| = number of documents
%
% * *|WO|* a 1 x |W| cell array of strings where |WO{k}| contains the kth
% vocabulary item and |W| is the number of distinct vocabulary items. Not
% needed for running the Gibbs sampler but becomes necessary when writing
% the resulting word-topic distributions to a file using the |writetopics|
% matlab function.
%
%%
% *OUTPUT*
%
% * *|WP|* a sparse matrix of size |W| x |T|, where |W| is the number of
% words in the vocabulary and |T| is the number of topics. |WP(i,j)|
% contains the number of times word |i| has been assigned to topic |j|.   
%
% * *|DP|* a sparse |D| x |T| matrix, where |D| is the number of documents.
% |DP(d,j)| contains the number of times a word token in document |d| has
% been assigned to topic |j|. 
%
% * *|Z|* a 1 x |N| vector containing the topic assignments where |N| is
% the number of word tokens. |Z(k)| contains the topic assignment for token
% |k|.    
%
%%
% *NOTES*
%
% If you have a text file with word-document counts, and would like to
% convert this text file into the matlab vectors |WS| and |DS|, use the
% conversion function |importworddoccounts|. This function assumes your text
% file is organized into three columns where each row contains the document
% index, the word index, and the word count. For example: 
% 1 2 10, 1 3 4, 2 2 6  (each comma representing a new line) 
% should be read as "word 2 occurs 10 times in doc 1, word 3 occurs 4 times
% in doc 1, and word 2 occurs 6 times in doc 2". (type help
% importworddoccounts for more information). To convert a text file of your
% vocabulary into a cell array of strings, use the "textread" function (a
% native Matlab function). For example, if your vocabulary is a text file
% "vocab.txt" with a different word on each row, then [ WO ] = textread(
% 'vocab.txt' , '%s' ) should convert this to an appropriate cell array of
% strings.
%
%%
%
%% AT (Author-Topic) MODEL
%
%%
% The input is a bag of word representation containing the number of times
% each words occurs in a document. Also needed is a matrix containing the
% authors present on each document. The outputs are the topic and author
% assignments to each word token as well as the counts of the number of
% times each word is assigned to a topic and the number of times a topic is
% assigned to an author.
%
%%
% *INPUT*
% 
% * *|WS|*  a 1 x |N| vector where |WS(k)| contains the vocabulary index of
% the kth word token, and |N| is the number of word tokens. The word
% indices are not zero based, i.e., min( |WS| )=1 and max( |WS| ) = |W| =
% number of distinct words in vocabulary
%
% * *|DS|*  a 1 x |N| vector where |DS(k)| contains the document index of
% the kth word token. The document indices are not zero based, i.e., min(
% |DS| )=1 and max( |DS| ) = |D| = number of documents
%
% * *|AD|*  a |A| x |D| sparse matrix where |A| is the number of distinct
% authors and |D| the number of documents. |AD(a,d)| = 1 when author |a| is
% present on document |d| and zero otherwise. 
%
% * *|WO|* a 1 x |W| cell array of strings where |WO{k}| contains the kth
% vocabulary item and |W| is the number of distinct vocabulary items. Not
% needed for running the Gibbs sampler but becomes necessary when writing
% the resulting word-topic distributions to a file using the |writetopics|
% matlab function.
%
% * *|AN|* a 1 x |A| cell array of strings where |AN{k}| contains the kth
% author name and |A| is the number of distinct authors. Not needed for
% running the Gibbs sampler but becomes necessary when writing the
% resulting author-topic distributions to a file using the |writetopics|
% matlab function.
%
%%
% *OUTPUT*
%
% * *|WP|* a sparse matrix of size |W| x |T|, where |W| is the number of
% words in the vocabulary and |T| is the number of topics. |WP(i,j)|
% contains the number of times word |i| has been assigned to topic |j|.   
%
% * *|AT|* a sparse |A| x |T| matrix, where |A| is the number of authors.
% |AT(a,j)| contains the number of times a word token associated with author |a| has
% been assigned to topic |j|. 
%
% * *|Z|* a 1 x |N| vector containing the topic assignments where |N| is
% the number of word tokens. |Z(k)| contains the topic assignment for token
% |k|.    
%
% * *|X|* a 1 x |N| vector containing the author assignments where |N| is
% the number of word tokens. |X(k)| contains the author assignment for token
% |k|. 
%
%%
%
%
%% LDA-COL (Collocation) MODEL
%
%%
% The input is a stream of words representation containing the order of
% words as they appear in documents. The outputs are the topic assignments
% to each word token and the assignments of tokens to a collocation or topic
% model. The output also includes the counts of the number of times each
% word is assigned to a topic and the number of times a topic is assigned
% to a document
%
%%
% *INPUT*
% 
% * *|WS|*  a 1 x |N| vector where |WS(k)| contains the vocabulary index of
% the kth word token, and |N| is the number of word tokens. The word
% indices are not zero based, i.e., min( |WS| )=1 and max( |WS| ) = |W| =
% number of distinct words in vocabulary. Note that the words are ordered
% according to occurence in documents, but that some words such as stop words can
% and should be omitted in this stream of words. 
%
% * *|DS|*  a 1 x |N| vector where |DS(k)| contains the document index of
% the kth word token. The document indices are not zero based, i.e., min(
% |DS| )=1 and max( |DS| ) = |D| = number of documents
%
% * *|WW|*  a |W| x |W| sparse matrix where |W(i,j)| contains the count of
% the number of times that word |i| follows word |j| in the word stream.
%
% * *|SI|*  a 1 x |N| vector where |SI(k)=1| only if the kth word can form
% a collocation with the (|k-1|)th word and |SI(k)=0| otherwise. Note that this
% representation is necessary because some consecutive words cross document
% boundaries and should not be allowed to form a collocation. Similarly,
% some words such as stop words in the original word stream might have been
% deleted. Any word at position |k| that follows a following a previously
% deleted word should have |SI(k)|=0  
%
% * *|WO|* a 1 x |W| cell array of strings where |WO{k}| contains the kth
% vocabulary item and |W| is the number of distinct vocabulary items. Not
% needed for running the Gibbs sampler but becomes necessary when writing
% the resulting word-topic distributions to a file using the |writetopics|
% matlab function.
%
%%
% *OUTPUT*
%
% * *|WP|* a sparse matrix of size |W| x |T|, where |W| is the number of
% words in the vocabulary and |T| is the number of topics. |WP(i,j)|
% contains the number of times word |i| has been assigned to topic |j|.   
%
% * *|DP|* a sparse |D| x |T| matrix, where |D| is the number of documents.
% |DP(d,j)| contains the number of times a word token in document |d| has
% been assigned to topic |j|. 
%
% * *|WC|* a 1 x |W| vector where |WC(k)| contains the number of times
% word |k| led to a collocation with the next word in the word stream.
%
% * *|C|* a 1 x |N| vector containing the topic/collocation assignments where |N| is
% the number of word tokens. |C(k)=0| when token |k| was assigned to the
% topic model. |C(k)=1| when token |k| was assigned to a collocation with word token |k|-1. 
%
% * *|Z|* a 1 x |N| vector containing the topic assignments where |N| is
% the number of word tokens. |Z(k)| contains the topic assignment for token
% |k|.    
%
%%
%
%% HMM-LDA MODEL
%
%%
% The input is a stream of words representation containing the order of
% words as they appear in documents. The outputs are the topic assignments
% to each word token and the assignments of tokens to a HMM state. The output 
% also includes the counts of the number of times each
% word is assigned to a topic and HMM state and the number of times a topic is assigned
% to a document. 
%
%%
% *INPUT*
% 
% * *|WS|*  a 1 x |N| vector where |WS(k)| contains the vocabulary index of
% the kth word token, and |N| is the number of word tokens. The word
% indices are not zero based, i.e., min( |WS| )=1 and max( |WS| ) = |W| =
% number of distinct words in vocabulary.  A word index of 0 denotes the
% end-of-sentence marker. Note that the words are ordered according to
% occurence in documents.
%
% * *|DS|*  a 1 x |N| vector where |DS(k)| contains the document index of
% the kth word token. The document indices are not zero based, i.e., min(
% |DS| )=1 and max( |DS| ) = |D| = number of documents
%
% * *|WO|* a 1 x |W| cell array of strings where |WO{k}| contains the kth
% vocabulary item and |W| is the number of distinct vocabulary items. Not
% needed for running the Gibbs sampler but becomes necessary when writing
% the resulting word-topic distributions to a file using the |writetopics|
% matlab function.
%
%%
% *OUTPUT*
%
% * *|WP|* a sparse matrix of size |W| x |T|, where |W| is the number of
% words in the vocabulary and |T| is the number of topics. |WP(i,j)|
% contains the number of times word |i| has been assigned to topic |j|.   
%
% * *|DP|* a sparse |D| x |T| matrix, where |D| is the number of documents.
% |DP(d,j)| contains the number of times a word token in document |d| has
% been assigned to topic |j|. 
%
% * *|MP|* a sparse |W| x |S| matrix where |S| is the number of HMM states. 
% |MP(i,j)| contains the number of times word |i| has been assigned to
% HMM state |j|. Note that HMM state 1 represents the LDA model and 2..|S|
% represent the "syntactic" HMM states
%
% * *|Z|* a 1 x |N| vector containing the topic assignments where |N| is
% the number of word tokens. |Z(k)| contains the topic assignment for token
% |k|. 
%
% * *|X|* a 1 x |N| vector containing the HMM state assignments where |N|
% is the number of word tokens. |X(k)| contains the assignment of the kth
% word token to a HMM state. Note that HMM state 1 represents the LDA model
% and 2..|S| represent the "syntactic" HMM states




