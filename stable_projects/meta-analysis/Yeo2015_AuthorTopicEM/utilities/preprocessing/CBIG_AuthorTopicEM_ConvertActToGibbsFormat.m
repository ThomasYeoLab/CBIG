function [WS, DS] = CBIG_AuthorTopicEM_ConvertActToGibbsFormat(act)
% [WS, DS] = CBIG_AuthorTopicEM_ConvertActToGibbsFormat(act)
%
% Convert activation data to appropriate format for the Gibbs sampler
%
% Input:
%   - act = array containing activation data across all experiments
% Output:
%   - WS and DS = arrays containing count required for the Gibbs sampler
%
% Example:
%   [WS, DS] = CBIG_AuthorTopicEM_ConvertActToGibbsFormat(act)
%   Produce appropriate input format for Gibbs sampler
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  num_words = sum(act, 2);
  WS = zeros(1, sum(act(:)));
  DS = zeros(1, sum(act(:)));
  
  count = 0;
  for i = 1:size(act, 1)
      DS(count+1:count+num_words(i)) = i;
      
      words_index = find(act(i, :) > 0);
      words_count = full(act(i, words_index));
      
      unique_word_count = unique(words_count);
      if(length(unique_word_count) > 1)
          error('CBIG_AuthorTopicEM_ConvertActToGibbsFormat assumes unique number of word count');
      end
      
      words = repmat(words_index, [unique_word_count 1]);
      words_count = sum(words_count);
      WS(count+1:count+words_count) = words(:)';
      
      count = count + num_words(i);
  end
