function [] = tirlVSexecTime(dateStart,dateEnd,connlocalDB)
% 
% Compares TIRL vs Execution time based on a date range
% excluding all members of Corsmed via their user_ids
% the function prints out the 2 different times in hours and seconds
% as well as the number of distinct users.
%
% INPUT
% dateStart = the starting date on format of 'yyyy-mm-dd'
% dateEnd   = the endind date on format of 'yyyy-mm-dd'
% connRemoteDB  = the connection object, necessary for connectio with DB
%
%



sqlquery=['SELECT exp.id,val.irl_info,exp.execTime_sec FROM edutool_globaldb.session_experiments as exp inner join edutool_globaldb.session_pulse_sequence as val on exp.id = val.exper_id and val.session_id = exp.session_id where val.dt_created > ''',dateStart,''' and val.dt_created <= ''',dateEnd,''' and exp.user_id != 1857 and exp.user_id!=790 and exp.user_id != 831 and exp.user_id != 927 and exp.user_id != 1139 and ', ...  
    'exp.user_id != 1857 and exp.user_id != 791 and exp.user_id !=2013 and exp.user_id != 933 and exp.user_id != 1929 and exp.user_id != 1438 and exp.user_id != 1543 and exp.user_id != 828 and exp.user_id != 2016 and val.irl_info IS NOT NULL and exp.execTime_sec IS NOT NULL and exp.status = ''finished'' ORDER BY val.irl_info ASC'];

%sqlquery=['SELECT val.irl_info, exp.execTime_sec FROM edutool_globaldb.session_experiments as exp inner join edutool_globaldb.session_pulse_sequence as val on exp.id = val.exper_id and val.session_id = exp.session_id where val.dt_created > ''',dateStart,''' and exps.user_id != 1857 and exps.user_id != 1438 and exps.user_id != 1543 and exps.user_id != 828 GROUP BY val.selected_value ORDER BY val.irl_info'];
sqlResults = exec(connlocalDB,sqlquery);
a = fetch(sqlResults);
results = a.Data;
totalTIRL = 0;
totalExecTime = 0;
for x=1:length(results)
    
    
    tirl = a.Data{x,2};
    if(contains(tirl,'m'))
        
        min = str2double(extractBefore(tirl,'m'));
        sec = extractBetween(tirl,':','s');
        sec = str2double(sec{1,1});
        
         totalTIRL = min*60+sec + totalTIRL;

    else
        totalTIRL =  str2double(extractBefore(tirl,'s')) + totalTIRL;
        
    end
    
    totalExecTime = str2double(a.Data{x,3}) + totalExecTime;
        
end


sqlquery2 = ['SELECT DISTINCT user_id FROM edutool_globaldb.session_experiments where dt_created > ''',dateStart,''' and dt_created <= ''',dateEnd,'''', ...
    'and user_id != 1857 and user_id!=790 and user_id != 831 and user_id != 927 and user_id != 1139 and user_id != 1857 and user_id != 933 and user_id != 1929 and user_id != 1438 and user_id != 1543 and user_id != 828 and execTime_sec IS NOT NULL and status = ''finished'''];
sqlResults2 = exec(connlocalDB,sqlquery2);
b = fetch(sqlResults2);
users = b.Data;



fprintf("Total Client Experiments : %d ",x)
fprintf("\nTotal TIRL time : %f seconds",totalTIRL);
fprintf("\nTotal TIRL time : %f hours \n",totalTIRL/3600);


fprintf("\nTotal Simulator time : %f seconds",totalExecTime);
fprintf("\nTotal Simulator time : %f hours \n",totalExecTime/3600);

fprintf("\nDifferent Users : %d \n",length(users));