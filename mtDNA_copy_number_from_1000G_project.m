sample_ids={'NA11831','NA10847','NA12872','NA12892','NA11931','NA11993','NA12812','NA12891','NA12761','NA12249','NA11995','NA11832','NA12776','NA12873','NA12045','NA12813','NA12155','NA12144','NA12716','NA11992','NA12003','NA11829','NA12006','NA10851','NA07037','NA12287','NA12154','NA12760','NA12414','NA12751','NA07346','NA12814','NA11881','NA12156','NA12750','NA12004','NA11830','NA11994','NA11840','NA12489','NA12234','NA12874','NA07347','NA07051','NA11920','NA11894','NA12815','NA12763','NA07000','NA12005','NA12762','NA12044','NA06985','NA06986','NA12043','NA12717','NA12828','NA11918','NA07357','NA06994'};

fileout='CEU_output.txt';


mt_gene.name={'MT-ATP6','MT-CO1','MT-CO2','MT-CO3','MT-CYB','MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6','MT-ATP8'};
mt_gene.pos=[8527,9207;5904,7445;7586,8269;9207,9990;14747,15887;3307,4262;4470,5511;10059,10404;10760,12137;10470,10766;12337,14148;14149,14673;8366,8572];

nc_gene.name={'AATF','ADAR','BSG','E2F4','MGAT1','NDUFA7','PGD','RBM8A','SOD1','TACC1','TALDO1','UBE2D2','YWHAB'};
nc_gene.chr={'17','1','19','16','5','5','1','1','21','8','11','5','20'};
nc_gene.pos=[35306175,35414171;154554533,154600475;571310,583493;67226068,67232821;180217541,180242652;180217541,180242652;10458649,10480233;145507598,145511444;33031935,33041244;38585704,38710546;747329,765024;138906016,139008019;43514317,43537173];

%%%%%%%%% prepare cmd line  %%%%%%%%%%%%
sample_ids=unique(sample_ids);

mt_pos=[];
for i=1:length(mt_gene.name)
    mt_pos=[mt_pos ' ' sprintf('MT:%d-%d',mt_gene.pos(i,:))];
end
mt_pos(1)='';

nc_pos=[];
for i=1:length(nc_gene.name)
    nc_pos=[nc_pos ' ' sprintf('%s:%d-%d',nc_gene.chr{i},nc_gene.pos(i,:))];
end
nc_pos(1)='';


%%
%%%%%%%%% get url_list for each individual
for i=1:length(sample_ids)
    url=sprintf('ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/%s/alignment/',sample_ids{i});
    try
        c=urlread(url);
    catch
        url_list{i}='';
        continue;
    end
    temp=regexp(c,'([^ \n]+)','tokens');
    if mod(length(temp),9) error('1'); end
    filename={};
    filesize=[];
    for j=0:length(temp)/9-1
        filename{j+1}=temp{j*9+9}{1};
        filesize(j+1)=str2num(temp{j*9+5}{1});
    end
    idx=find(filesize==max(filesize));
    if ~strcmp(filename{idx}(end-3:end),'.bam') error('2'); end
    if ~strcmp(filename{idx+1},sprintf('%s.bai',filename{idx})) error('3'); end
    url_list{i}=sprintf('%s%s',url,filename{idx});
    url_2=sprintf('%s%s',url,filename{idx+1});
    %urlwrite(url_2,filename{idx+1});
end

%%%%%%%% output
%%
fid=fopen(fileout,'w');
fprintf(fid,'nc_gene');
for i=1:length(sample_ids)
    fprintf(fid,'\t%s',sample_ids{i});
end
fprintf(fid,'\n');

fprintf(fid,'MT');
for i=1:length(sample_ids)
    if isempty(url_list{i})
        fprintf(fid,'\tNaN');
        continue;
    end
    disp(sprintf('MT: sample # %d',i));
    cmd=sprintf('samtools view %s %s|wc -l',url_list{i},mt_pos);
    [a l]=system(cmd);
    l(end)='';
    
        fprintf(fid,'\t%s',l);
    
end
fprintf(fid,'\n');

fprintf(fid,'nc_gene');
for i=1:length(sample_ids)
    if isempty(url_list{i})
        fprintf(fid,'\tNaN');
        continue;
    end
    disp(sprintf('Gene: sample # %d',i));
    cmd=sprintf('samtools view %s %s|wc -l',url_list{i},nc_pos);
    [a l]=system(cmd);
    l(end)='';
    
    
        fprintf(fid,'\t%s',l);
    
end
fprintf(fid,'\n');

fclose(fid);