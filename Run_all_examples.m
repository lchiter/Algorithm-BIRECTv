clear all;
global  example
fid = fopen('results.txt','wt');
fevaluation=[];
   fprintf(fid,'%s\t','No');
    fprintf(fid,'%s\t','f_star');
    fprintf(fid,'%s\t\t','f(xmin)');
    fprintf(fid,'%s\t','f_eval');
    fprintf(fid,'%s\t','it');
    fprintf(fid,'%s\t','CPU');   
    fprintf(fid,'\n');

for example = 1:1%[1:10 15:36 38 ]
    Birect2
    fevaluation=[fevaluation ;example fmin f_eval];
    fprintf(fid,'%d\t',example);
    fprintf(fid,'%15.8e\t',f_star);
    fprintf(fid,'%15.8e\t',fmin);
    fprintf(fid,'%d\t',f_eval);
    fprintf(fid,'%d\t',it);
%    fprintf(fid,'%d\t',CPU);
    fprintf(fid,'\n');
end
fevaluation;
fprintf(fid,'\n');
fprintf(fid,'%s\t\t   ','Average');fprintf(fid,'%10.3f',mean(fevaluation(:,3)));
fprintf(fid,'\n');
fprintf(fid,'%s\t\t   ','Median');fprintf(fid,'%10.3f',median(fevaluation(:,3)));

fclose(fid);