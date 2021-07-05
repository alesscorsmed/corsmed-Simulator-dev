function [R, S] = redisSetJson(R, key, value)

S = 'OK';

if ~isstr(value)
  S = 'ERROR - SET VALUE MUST BE A STRING';
  return
end

if ~strcmp(R.status, 'open')
  S = 'ERROR - NO CONNECTION';
  return
end

[Response, R, S] = tools.redis.redisCommand(R, tools.redis.redisCommandJson(sprintf('SET %s', key),value));

if strcmp(Response, '+OK')
  S = Response;
  return
end

if contains(Response,'WrongBufferSize')
    S = Response;
    return
end