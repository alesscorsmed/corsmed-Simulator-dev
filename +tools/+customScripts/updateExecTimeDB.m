function updateExecTimeDB(dateStart,connlocalDB)
%
% Function used to update the exec_time of experiments in local DB
% based on JSON files (STATS). Reads the JSON Files from the EFS and finds 
% corresponding experiments based on user,pulse,experiment and date IDs.
%
% INPUT 
% dateStart = the datetime from which the update of the exec_time takes
%             place. Format used is 'yyyy-mm-dd'
%
% connlocalDB = the connection object, used for connecting to the local DB.
%
sqlquery = ['SELECT DISTINCT user_id FROM edutool_globaldb.session_experiments where dt_created > ''',dateStart,''' and status = ''finished''  '];
sqlResults = exec(connlocalDB,sqlquery);
a = fetch(sqlResults);
users = a.Data;

if(~strcmp(users{1},'No Data'))
    for x=1:length(users)

        sqlquery = ['SELECT id,pulseq_id,dt_created FROM edutool_globaldb.session_experiments where user_id=''',num2str(users{x}),''' and dt_created > ''',dateStart,''' and status=''finished'' '];
        sqlResults = exec(connlocalDB,sqlquery);
        b = fetch(sqlResults);
        experiments = b.Data;

        for i=1:size(experiments,1)

        experstr = "E"+experiments{i,1}+"_";
        pulsestr = "P"+experiments{i,2}+"_";
        
        dt_created = experiments{i,3};
        match = [" ",":","-","."];
        datetimestr = erase(dt_created,match);
        datetimestr = datetimestr(1:end-3);
        userstr = "U"+users(x);

        P = '/efs-mount-point/S20/STATS/edutool';
        S = dir(fullfile(P,'*.json'));

            for k = 1:numel(S)
                N = S(k).name;
                R1 = regexp(N,userstr,'match','all'); 
                R2 = regexp(N,experstr,'match','all'); 
                R3 = regexp(N,pulsestr,'match','all');
                R4 = regexp(N,datetimestr,'match','all');
                if(~isempty(R2) && ~isempty(R1) && ~isempty(R3) && ~isempty(R4))
                    jsonFile = S(k);

                        fullfilepath = fullfile(jsonFile.folder,jsonFile.name);
                        fid = fopen(fullfilepath);
                        raw = fread(fid,inf);
                        str = char(raw');
                        fclose(fid);
                        data = jsondecode(str);
                        disp(data.timeTotal)

                        sqlquery = ['UPDATE edutool_globaldb.session_experiments SET execTime_sec = ''',num2str(data.timeTotal),''' where user_id=''',num2str(users{x}),''' and id = ''',num2str(experiments{i,1}),''' and pulseq_id = ''',num2str(experiments{i,2}),''' and dt_created > "2021-03-20"'];
                        sqlResults = exec(connlocalDB,sqlquery);

                end
            end


        end

    end
else
    disp("No Experiments with NULL execution time were found.")  
end
