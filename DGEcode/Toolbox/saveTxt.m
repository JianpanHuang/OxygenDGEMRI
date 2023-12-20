function saveTxt(pathStr,saveMat)
% save 'saveMat' as a txt file to path 'pathStr'
fid=fopen(pathStr,'wt');
smat = saveMat;
[m,n]=size(smat);
 for i=1:1:m
   for j=1:1:n
      if j==n
        fprintf(fid,'%g\n',smat(i,j));
     else
       fprintf(fid,'%g\t',smat(i,j));
      end
   end
end
fclose(fid);