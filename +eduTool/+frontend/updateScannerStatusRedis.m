function [c, status] = updateScannerStatusRedis(redisConn,message,key)

keyNew      = [key,'_status'];
[c, status] = tools.redis.redisSet(redisConn,keyNew,message);