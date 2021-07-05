function R = redisDisconnect(R)

fclose(R);

disp('Closing Redis connection')
