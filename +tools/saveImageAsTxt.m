function saveImageAsTxt(normalizedImageResized,txt_name_image)
% Save complex image as a .txt file
fileID = fopen(txt_name_image,'w');
if isreal(normalizedImageResized)
    for imageLine=1:size(normalizedImageResized,1)
        realk = normalizedImageResized(imageLine,:);
        realk = realk(:);
        fprintf(fileID,'%f;',realk(1:end-1).');
        fprintf(fileID,'%fi|',realk(end).');
    end
else    
    for imageLine=1:size(normalizedImageResized,1)
        realk = real(normalizedImageResized(imageLine,:));
        imagk = imag(normalizedImageResized(imageLine,:));
        realk = realk(:);
        imagk = imagk(:);
        fprintf(fileID,'%f%+fi;',[realk(1:end-1),imagk(1:end-1)].');
        fprintf(fileID,'%f%+fi|',[realk(end),imagk(end)].');
    end
end
fclose(fileID);