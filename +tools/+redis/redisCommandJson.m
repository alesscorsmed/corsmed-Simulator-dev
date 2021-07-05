function cmd = redisCommandJson(Cmd,jsonStr)

Cmd = strtrim(Cmd);
words = regexp(Cmd, '\s', 'split');
words{end+1} = jsonStr;

n_words = numel(words);

cmd = sprintf('*%d', n_words);
for ix = 1 : n_words
  word = words{ix};
  word_length = numel(word);
  cmd = sprintf('%s\n$%d\n%s', cmd, word_length, word);
end

% note the final \r\n is added by the port
