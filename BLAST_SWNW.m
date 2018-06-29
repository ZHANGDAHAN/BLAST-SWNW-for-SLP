function [output]=BLAST_SWNW(varargin)
%   The funtion is to carrying out BLAST+SWNW
%   function [output]=BLAST_SWNW(input_file,drill,output_pwd,Blastpath) can
%   give user the best subcellular prediction result of inputfile by
%   BLAST+SWNW. In output, there are three vetor：Protein_Id, Prediction
%   result ,and P value.
%   Examples:
%
%       Return output by inputting 4 variables
%       [output] = BLAST_SWNW(input_file,drill,output_pwd,Blastpath)
tic;
narginchk(4,inf)
if length(varargin)==4
    %% save innital path
    input_file=varargin{1};
    drill=varargin{2};
    output_pwd=varargin{3};
    Blastpath=varargin{4};
    %     innital_path=pwd;
    cd (Blastpath)
    %% construct a drill
    string_need_run=['system(''makeblastdb.exe -in ',drill,' -dbtype prot -out ',drill,''')'];
    eval(string_need_run)
    %% run BLAST
    outputfile=strcat(output_pwd,'\output.blast');
    string_need_run=['system(''blastp.exe -query ',input_file,' -out ',outputfile,'  -db ',drill,'  -num_descriptions 10  -num_threads 8 -outfmt 6 -evalue 10 '')'];
    eval(string_need_run)
    %% read blastp output
    string_need_run=['data=blastreadlocal(''',outputfile,''',8);'];
    eval(string_need_run)
    
    %% run 'SWNW'
    testseqnum = length( data );
    data_Query = cell(testseqnum,1);
    for i=1:testseqnum
        data_Query{i,1}=data(i).Query;    %提取出data中的Query
        %         location_real(1,i)=str2num(data_Query{i,1}(strfind(data_Query{i,1},'#')+6:strfind(data_Query{i,1},'_')-1));
    end
    
    % align two sequences using Smith-Waterman algorithm and Needleman-Wunsch algorithm
    string_need_run=strcat('[ENTRY, seq]=fastaread(''',input_file,''');');
    [Header_drill, seq_drill]=fastaread(drill);
    eval(string_need_run);
    [ min_P, ~, best_loc ] = sw_nw_dahan( ENTRY,  seq,  data, Header_drill, seq_drill );
    output=cell(testseqnum,3);
    %% predict output
    classify_24={'Cellmembrane', 'cellwall','chloroplast','chloroplast membrane','chloroplast thylakoid','chloroplast thylakoid membrane','Cytoplasm','Endoplasmic reticulum','Endoplasmicreticulummembrane','Golgiapparatus','Golgistackmembrane','Lysosome','Lysosome membrane','Mitochondrion','Mitochondrion membrane','Nucleus','Nucleus membrane','Peroxisome','Peroxisome membrane','Plastid','Plastid membrane','Secreted','Vacuole','Vacuole membrane'};
    for i=1:testseqnum
        output{i,1}=data(i).Query;
        output{i,2}=classify_24{str2num(best_loc{i}((strfind(best_loc{i},'class')+5):(strfind(best_loc{i},'_')-1)))};
        output{i,3}=min_P(i);
        
    end
    %     run_time=toc;
    % else   %%% return Error message.
    %     error('Error occurred: no enough input! \nprogram failed.',class(7))
    %
end

end

% dataset4_20chong_blastp_swnw.m 特别定制版
function [ min_P,loc_min_P,best_loc ] = sw_nw_dahan( ENTRY,  seq,  data_filter_1en1, Entry_drill, seq_drill )
%sw_nw_dahan 用于对blast_去冗余后的结果进行swnw处理，返回最佳定位和最优结果和输入序列的P值（P=p_nw*p_sw;）
%   此处显示详细说明
min_P=zeros(length(data_filter_1en1),1);
loc_min_P=zeros(length(data_filter_1en1),1);
best_loc=cell(length(data_filter_1en1),1);
k=0.20441;
lanmuda=0.16931;
for i=1:length(data_filter_1en1)
    %     i
    id_origin=data_filter_1en1(i).Query;
    w1=find(ismember(ENTRY,id_origin)>0);
    seq1=seq{w1};
    seq1(isspace(seq1)) = [];%去掉空格,只保留20种氨基酸的字符
    seq1=strrep(seq1,'B','');  % B，U，J，O,X，Z
    seq1=strrep(seq1,'J','');
    seq1=strrep(seq1,'O','');
    seq1=strrep(seq1,'U','');
    seq1=strrep(seq1,'X','');
    seq1=strrep(seq1,'Z','');
    p_sw=zeros(1,1);p_nw=zeros(1,1);P=zeros(1,1);
    for j=1:length(data_filter_1en1(i).Hits)
        %返回第j条序列的
        id=data_filter_1en1(i).Hits(j).Name;
        %        t=strfind(id,'#');
        entry_predict=id; %(1:t(1)-1);            %%%%提取blastp结果中的entry name
        w2=find(ismember(Entry_drill,entry_predict)>0);
        seq2=seq_drill{w2};
        seq2(isspace(seq2)) = [];  %去掉空格  'ACDEFGHIKLMNPQRSTVWY' B，U，J，O,X，Z
        seq2=strrep(seq2,'B','');  % B，U，J，O,X，Z
        seq2=strrep(seq2,'J','');
        seq2=strrep(seq2,'O','');
        seq2=strrep(seq2,'U','');
        seq2=strrep(seq2,'X','');
        seq2=strrep(seq2,'Z','');
        [Score_sw, ~] = swalign(seq1, seq2);
        [Score_nw, ~] = nwalign(seq1, seq2);
        m=length(seq1);n=length(seq2);
        p_sw(1,j)=1-exp(-k*m*n*exp(-lanmuda*Score_sw));
        p_nw(1,j)=1-exp(-k*m*n*exp(-lanmuda*Score_nw));
        P(1,j)=p_nw(1,j)*p_sw(1,j);
    end
    if length(data_filter_1en1(i).Hits)>0
        [min_P(i,1),loc_min_P(i,1)]=min(P);
        %     best_loc{i,1}=data_filter_1en1(i).Query;
        best_loc{i,1}=data_filter_1en1(i).Hits(loc_min_P(i,1)).Name;
    else
        min_P(i,1)=nan;
        loc_min_P(i,1)=nan;
        
        
        
        
        
        
        best_loc{i,1}='';
    end
end
end











