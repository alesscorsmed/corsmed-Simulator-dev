%% Connect to redis
try
    ctrl = mps.cache.control('myRedisConnection','Redis',...
    'Host','redis-nodeport.internal.integration.corsmed.com',...
    'Port',30007);
    start(ctrl)
catch ME
    
end
c = mps.cache.connect('myCache','Connection','myRedisConnection');

disp('Redis connection established.')

%% Retrieve data from redis
basicText = 'anatomicalMIDAspatial_';
spatial = zeros(80640000,3);

batchSize = 1000000;
noOfFiles = 81;
tic
for i=1:noOfFiles
    startIndex = (i-1)*batchSize+1;
    if i == noOfFiles
        endIndex = size(spatial,1);
    else
        endIndex = i*batchSize;
    end
    
    key = [basicText,num2str(i)];
    spatial(startIndex:endIndex,:) = get(c,key);
    disp(['Part ',num2str(i),' retrieved successfully'])
end
toc