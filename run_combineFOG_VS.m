% Voor samenvoegen van annotaties en corrigeren
clear all; close all;

Ann_folder1='F:\Vibrating Socks FOG Scoring\Annotatie - Emilie\';
Ann_folder3='F:\Vibrating Socks FOG Scoring\Annotatie - Marleen\';
folder_combined='F:\Vibrating Socks FOG Scoring\Output ElanPaper\';

subjects1=dir(fullfile(Ann_folder1, '*VS*'));

rater=3;    % 2 jamie, 3 marleen

% Subjects aanpassen voor rater 2 of 3
if rater == 2
    %subjects=dir(fullfile(Ann_folder2, '*VS*'));
    Ann_folderrater2=Ann_folder2;
elseif rater ==3
    Ann_folderrater2=Ann_folder3;
end

subjects=dir(fullfile(Ann_folderrater2, '*VS*'));
%%
for VSnmmr2=8
    
    [Video_files1, Video_files2] = createVidlist(VSnmmr2, Ann_folder1,...
        subjects1, rater ,Ann_folderrater2, subjects);
    indxVid = linkvids(Video_files1, Video_files2);
    %%
    for filenummer1=17 %1:length(Video_files1)
        
        filename_rater1=[Video_files1(filenummer1).folder,'\',Video_files1(filenummer1).name];
        filenummer2=find (indxVid == filenummer1);
        %%
        % Niet alle files zijn geannoteerd door Ann2, wanneer duidelijk geen freezing voor Ann1
        % Deze worden dus overgeslagen
        if ~isempty(filenummer2)
            filename_rater2=[Video_files2(filenummer2).folder ,'\',Video_files2(filenummer2).name];
            
            bestandsnaam=Video_files1(filenummer1).name;
            bestandsnaam=bestandsnaam(1:end-7);
            
        filename_combined=fullfile(folder_combined, bestandsnaam);
        agreement_table=fullfile(folder_combined, 'agreement_table.tsv');
        
        combine_FOGannotations(filename_rater1, filename_rater2, filename_combined,...
            agreement_table, subjects1(VSnmmr2).name, 'include', 2)

        end       

    end
end
%%
agreement_Tabel=readtable(agreement_table, 'FileType', 'text', 'ReadVariableNames', 1, 'HeaderLines', 0);
%% Lokale functies

function [Video_files1, Video_files2] = createVidlist(VSnmmr2, Ann_folder1, subjects1, rater ,Ann_folderrater2, subjects)
%Voor debuggen
% VSnmmr2=...;

%%
VSname={subjects(VSnmmr2).name};
indxsub1 = strcmpi ({subjects1.name}, VSname);  % Vind welke overeenkomt met map Emilie
VSnmmr1 = find(indxsub1);

% Map van Emilie per meetmoment
measurementdates=dir(fullfile([Ann_folder1, subjects1(VSnmmr1).name],...
    '*20*')); % Map heeft altijd 20 (2020 of 2021) er in.

filespath1a=[Ann_folder1, subjects1(VSnmmr1).name, '\' ,measurementdates(1).name, '\'];
Video_files1a=dir(fullfile(filespath1a, '*.txt*'));

if length(measurementdates)==2
    filespath1b=[Ann_folder1, subjects1(VSnmmr1).name, '\' ,measurementdates(2).name, '\'];
    Video_files1b=dir(fullfile(filespath1b, '*.txt*'));
    Video_files1=[Video_files1a; Video_files1b];
elseif  length(measurementdates) ==1
    Video_files1=Video_files1a;
end


if rater == 2           % Map Jamie
    Video_files2=dir(fullfile([Ann_folderrater2, subjects(VSnmmr2).name],...
        '*.txt*'));
elseif rater == 3       % Map Marleen
    filespath2a=[Ann_folderrater2, subjects(VSnmmr2).name, '\' ,measurementdates(1).name, '\'];
    Video_files2a=dir(fullfile(filespath2a, '*.txt*'));
    if length(measurementdates)==2
        filespath2b=[Ann_folderrater2, subjects(VSnmmr2).name, '\' ,measurementdates(2).name, '\'];
        Video_files2b=dir(fullfile(filespath2b, '*.txt*'));
        Video_files2=[Video_files2a; Video_files2b];
    elseif  length(measurementdates) ==1
        Video_files2=Video_files2a;
    end
    
end
end

% Linkvids; link 2 lijsten met videos van ongelijke aantallen aan elkaar
function indxVid = linkvids(Video_files1, Video_files2)
list1=char(Video_files1.name);
list2=char(Video_files2.name);

for M=1:length(Video_files1)
    vidname1=strsplit(list1(M,:), '_');
    vidname1=strrep(vidname1(2),' ','');
    videos1{M}=vidname1;
end

for L=1:length(Video_files2)
    vidname2=strsplit(list2(L,:), '_');
    videos2{L}=vidname2(2);
end

% Alleen video nummers
videos1=str2double(string(videos1));
videos2=str2double(string(videos2));

indxVid=0; % pre-allocate, ook voor wanneer geen enkele overeen komt

% Check wanneer deze overeen komen
for K=1:length(Video_files2)
    test1 = find (videos1 == videos2(K));
    if ~isempty(test1)      % als test1 niet leeg is, dan er naar toe schrijven
        indxVid(K)=test1;
    end
    clear test1
end
end