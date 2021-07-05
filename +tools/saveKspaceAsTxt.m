function saveKspaceAsTxt(normalizedImageResized,txt_name_kSpace,FT)
% Save kspace as a .txt file for the kspace tool
normalizedImageResized_kspace = FT(normalizedImageResized);
fileID = fopen(txt_name_kSpace,'w');
for kspaceLine=1:size(normalizedImageResized_kspace,1)
    realk = real(normalizedImageResized_kspace(kspaceLine,:));
    imagk = imag(normalizedImageResized_kspace(kspaceLine,:));
    realk = realk(:);
    imagk = imagk(:);
    fprintf(fileID,'%f%+fi;',[realk(1:end-1),imagk(1:end-1)].');
    fprintf(fileID,'%f%+fi|',[realk(end),imagk(end)].');
end
fclose(fileID);