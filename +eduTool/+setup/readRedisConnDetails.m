function redisConnDetails = readRedisConnDetails(txtInput)

fid = fopen(txtInput,'r');
urlstring = textscan(fid,'%s');
fclose(fid);

IndexRedisAddress   = strfind(urlstring{1}, '#address');
isOne               = cellfun(@(x)isequal(x,1),IndexRedisAddress);
row_RedisAddress    = find(isOne);

IndexRedisPort      = strfind(urlstring{1}, '#port');
isOne               = cellfun(@(x)isequal(x,1),IndexRedisPort);
row_RedisPort       = find(isOne);

IndexInputBuffer    = strfind(urlstring{1}, '#inputBufferSize');
isOne               = cellfun(@(x)isequal(x,1),IndexInputBuffer);
row_InputBuffer     = find(isOne);

IndexOutputBuffer   = strfind(urlstring{1}, '#outputBufferSize');
isOne               = cellfun(@(x)isequal(x,1),IndexOutputBuffer);
row_OutputBuffer    = find(isOne);


redisConnDetails.address            = char(urlstring{1}(row_RedisAddress+1,1));
redisConnDetails.port               = char(urlstring{1}(row_RedisPort+1,1));
redisConnDetails.inputBufferSize    = char(urlstring{1}(row_InputBuffer+1,1));
redisConnDetails.outputBufferSize   = char(urlstring{1}(row_OutputBuffer+1,1));