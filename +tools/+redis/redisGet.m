function [Value, R, S] = redisGet(R, key)

S = 'OK';
Value = [];

if ~strcmp(R.Status, 'open')
  S = 'ERROR - NO CONNECTION';
  return;
end

[Response, R, S] = tools.redis.redisCommand(R, tools.redis.redisCommandString(sprintf('GET %s', key)));

if Response(1) == '-'
  S = Response;
  return
end

% response $-1 means nonexistant key
if Response(2) == '-'
  S = 'ERROR - NONEXISTANT KEY';
  return
end

if contains(Response,'WrongBufferSize')
    S = Response;
    return
end

if Response(1) ~= '$'
  S = Response;
  return
end

Value = tools.redis.redisParseBulkReply(Response);