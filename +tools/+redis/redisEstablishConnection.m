function R = redisEstablishConnection(redisAddress,redisPort,...
    inputBufferSize,outputBufferSize)
% This connection is utilized using the redis client for MATLAB (found
% here: github.com/dantswain/redis-matlab

if nargin < 3
    inputBufferSize     = 10000000; % 10MB
    outputBufferSize    = 10000000; % 10MB 
end

[R] = tools.redis.redisConnection(redisAddress,str2num(redisPort),...
    inputBufferSize,outputBufferSize);

% Ping redis server to see if connection is established
[Response,~,~] = tools.redis.redisPing(R);
disp(['Response from redis: ',Response])
% Redis vill return 'PONG' if success
if ~strcmp(Response, 'PONG')
    tools.redis.redisDisconnect(R)
    ME = MException('eduTool:failedRedisConnection',...
        '%s : failed Redis connection',functionName);
    throw(ME);
else
    disp('Connection with Redis just established')
end