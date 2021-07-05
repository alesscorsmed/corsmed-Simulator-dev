% This is the way how we can load environment variables instead of reading
% from files
redisAddress    =  'redis-nodeport.internal.integration.corsmed.com'; % getenv('CORSMED_REDIS_ADDRESS');
redisPort       =	30007; % str2double(getenv('CORSMED_REDIS_PORT'))
% redisKey     =	getenv('CORSMED_REDIS_SECRET_KEY');

% Connect to a redis server
[R] = tools.redis.redisConnection(redisAddress, redisPort);

% Ping redis server to see if connection is established
[Response, C2, Status2] = tools.redis.redisPing(R);

% Redis vill return '+OK' if success and if does, for this purpose try to
% load model from the cache
if strcmp(Response, 'PONG')
    key = 'COURSE_3_790_6';
%     key = 'COURSE_ID_147.mat';
%     key = '20210118_test';
%     value = '3';
    
    % Add file to a redis cache 
%     [C3, Status3] = tools.redis.redisSet(R, key, value);
    
    % Retrieve file from redis
    [model,R,S] = tools.redis.redisGet(R, key);
end

% Disconnect from the server
tools.redis.redisDisconnect(R)