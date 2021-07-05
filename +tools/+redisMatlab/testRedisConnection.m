try
    ctrl = mps.cache.control('myRedisConnection','Redis',...
        'Host','redis-nodeport.internal.integration.corsmed.com',...
        'Port',30007);
    start(ctrl)
catch ME
    
end

c = mps.cache.connect('myCache','Connection','myRedisConnection');
put(c,'count',33);
isKey(c,'count')
x = get(c,'count');

disp(num2str(x))

%%
% % get all keys from cache
keys(c)
% 
% % remove key from cache
% remove(c,'count2')

