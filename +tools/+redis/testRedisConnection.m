redis.address   = 'cache.internal.integration.corsmed.com';
redis.port      = '6379';
keyInput        = 'EXPERIMENT_1_1';
keyOutput       = 'CX_20200201';

inputBufferSize     = 10000000;
outputBufferSize    = 10000000;

% Establish connection
R = tools.redis.redisEstablishConnection(redis.address,redis.port,...
    inputBufferSize,outputBufferSize);

%% GET
[jsonStructure,R] = tools.redis.redisGetJsonWrapper(R,keyInput);
% tools.redis.redisDisconnect(R);

disp(' ')

%% SET
R = tools.redis.redisSetJsonWrapper(R,keyOutput,jsonStructure);

%% GET ALL KEYS
[Response, R, S] = tools.redis.redisCommand(R,...
    tools.redis.redisCommandString(sprintf('keys *')));
tools.redis.redisDisconnect(R);

ResponseCell    = regexp(Response, '\s', 'split');
redisKeys       = sort(ResponseCell(1,5:4:end));

fprintf(1, '\n  REDIS KEYS : \n');
for i=1:size(redisKeys,2)
    fprintf(1,'%s\n',redisKeys{1,i});
end

%% Close connection
tools.redis.redisDisconnect(R);