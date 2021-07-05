function [] = batchReadFile(filename,tab_name,courseID,connLocalDB)
%
% TOOLS.BATCHREADFILE
%
%	Runs a experiment.
%
% INPUT
%   filename            excel file path 
%   tab_name            the excel's file tab 
%   courseID           	id of the course
%   connLocalDB         local connection object
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%

%Pulse Parameter Rows must start from row number 4.
[num,txt,raw] = xlsread(filename,tab_name);

for i=1:length(raw(:,1))
    if(~isnan(raw{i,1}))
        index = i;
        break;
    end
end

num_exp = length(txt(:,1))-2;

exec(connLocalDB, ...
['SET @out_batch_id = 0;']);

 exec(connLocalDB,...
['CALL edt_tool_local.load_batch(''',num2str(courseID),''',@out_batch_id)']);

%         fetchedbatch = exec(conn_localdb,...
%         ['SELECT @out_batch_id']);
%         batchout = fetch(fetchedbatch);

exec(connLocalDB, ...
['INSERT INTO edt_tool_local.batch_load (course_id,file_id)',...
'VALUES (''',num2str(courseID),''',''',filename,''')']);

batchload_id = exec(connLocalDB,...
['SELECT MAX(id) FROM batch_load']);
batchload_id = fetch(batchload_id);
batchload_id = batchload_id.Data{1};


for i=1:num_exp
    excel_row = raw(i+2,2:5);
    disp(excel_row);
    
    pulseq_name = extractBefore(raw(i+2,1),'_');
    pulseq_name = pulseq_name{1};
    pulseq_params = excel_row{1,1};
    
    points_params = excel_row{1,2};
    simpoints_params = excel_row{1,3};
    attribute_params = excel_row{1,4};
    
    exec(connLocalDB, ...
            ['INSERT INTO edt_tool_local.batch_protocol (pulseq_name,pulseq_params,spatial_points,spatial_sim_points,configurations,batch_id)',...
            ' VALUES (''',pulseq_name,''',''',pulseq_params,''',''',points_params,''',''',simpoints_params,''',''',attribute_params,''',''',num2str(batchload_id),''')']);

            pulseq_params = '';
            simpoints_params = '';
            points_params = '';
            attribute_params = '';

end

try
    
    exec(connLocalDB, ...
            ['SET @out_arg = 0;']);

      exec(connLocalDB,...
    ['CALL edt_tool_local.batch_run(''',num2str(batchload_id),''',@out_arg)']);

catch ME    
    parameterString = "Failed to Load Experiment";
    disp(parameterString);
end

parameterString = "Success Loading Batch";
disp(parameterString);