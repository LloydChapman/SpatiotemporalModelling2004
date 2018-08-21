function PrintModesAndCIsToFile(x,str)

fid = fopen(str,'w');
for i=1:size(x,1)
    fprintf(fid, ['& %.2f & (%.2f--%.2f) & %.0f & (%.0f--%.0f) ',repmat('& %.1e & (%.1e--%.1e) ',1,2),'& %.1f \\\\ \n'], x(i,:)');
end
fclose(fid);