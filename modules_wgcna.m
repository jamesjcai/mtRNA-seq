clear;
clc;
fclose all;
load ..\data_gene.mat;
z=data_gene(:,2);
filename='modules.txt';
fid=fopen(filename);
content=textscan(fid,'%s%s','delimiter','\t','headerLines',1);
gene=content{1};
module=content{2};
fclose(fid);

[module_uni,ia,ic]=unique(module);

fout=fopen('modules_table.txt','w');
fprintf(fout,'%s\t%s\n','modules','gene');
 for i=1:length(module_uni)
     idx=find((ic==i));
     fprintf(fout,'%d\t',i);
     for j=1:length(idx)-1
         y=gene{idx(j)};
         idy=find(ismember(data_gene(:,1),y));          
         if sum(ismember(data_gene(:,1),y))>0;
             y=z{(idy)};
         end
         
     fprintf(fout,'%s ,',y);
     
     end
      fprintf(fout,'%s\n',gene{idx(j)+1});
 end
