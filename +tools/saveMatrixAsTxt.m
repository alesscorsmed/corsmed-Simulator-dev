function saveMatrixAsTxt(matrix,txtSaveFileName)
% Save complex image as a .txt file
fileID = fopen(txtSaveFileName,'w');
if isreal(matrix)
    for imageLine=1:size(matrix,1)
        realk = matrix(imageLine,:);
        realk = realk(:);
        fprintf(fileID,'%f;',realk(1:end-1).');
        fprintf(fileID,'%fi|',realk(end).');
    end
else    
    for imageLine=1:size(matrix,1)
        realk = real(matrix(imageLine,:));
        imagk = imag(matrix(imageLine,:));
        realk = realk(:);
        imagk = imagk(:);
        fprintf(fileID,'%f%+fi;',[realk(1:end-1),imagk(1:end-1)].');
        fprintf(fileID,'%f%+fi|',[realk(end),imagk(end)].');
    end
end
fclose(fileID);