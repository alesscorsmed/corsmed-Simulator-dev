function c = redisEstablishConnectionMatlab(redisAddress,redisPort)

try
    ctrl = mps.cache.control('myRedisConnection','Redis',...
        'Host',redisAddress,...
        'Port',str2num(redisPort));
    start(ctrl)
catch ME
    disp('Redis connection already exists')
end

% c = mps.cache.connect('myCache','Connection','myRedisConnection')
c = mps.cache.connect('myCache')