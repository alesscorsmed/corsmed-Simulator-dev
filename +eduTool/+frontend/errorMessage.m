function edt_db_errorMessage(conn,pulseq_id,error_message)

curs2 = exec(conn, ...
    ['UPDATE edt_tool_local.pulse_sequence SET info=''',error_message,...
    ''' WHERE id=',num2str(pulseq_id)]);