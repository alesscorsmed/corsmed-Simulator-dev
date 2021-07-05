function [R, S] = redisDEL(R, key)

if ~strcmp(R.status, 'open')
  S = 'ERROR - NO CONNECTION';
  return
end

[Response, R, S] = tools.redis.redisCommand(R, tools.redis.redisCommandString(sprintf('DEL %s', key)));

if strcmp(Response(1),':') && strcmp(Response(2),'1')
    S = 'OK';
end
