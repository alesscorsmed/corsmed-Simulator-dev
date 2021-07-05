function updateScannerStatus(expControl,message,type,progress)

if nargin==1
    message     = 'Experiment is running';
    type        = 'update';
    progress    = '';
elseif nargin==2
    type        = 'update';
    progress    = '';
elseif nargin==3
    progress = '';
end

if nargin<4 && isfield(expControl,'progress')
    progress = num2str(expControl.progress);
end

if isfield(expControl,'connLocalDB')
    conn = expControl.connLocalDB;

    if ~isempty(conn)
        exec(conn, ...
            ['UPDATE edt_tool_local.global_configuration SET selected_value=''',message,...
            ''' WHERE name=''run_info''']);
    end
else
    tools.updateJsonScannerProgress(expControl,type,message,progress);
end