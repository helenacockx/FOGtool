%% set-up:
clear all; close all;

folder_rater1='\\dcn-srv.science.ru.nl\dcn\biophysics\prompt\freezing_fnirs\data\processed\annotations\Helena';
folder_rater2='\\dcn-srv.science.ru.nl\dcn\biophysics\prompt\freezing_fnirs\data\processed\annotations\Yuli';
folder_combined='\\dcn-srv.science.ru.nl\dcn\biophysics\prompt\freezing_fnirs\data\processed\annotations\combined';

subjects={'PD10', 'PD11'};

sf=1000; % choose sampling frequency
ts=(1/sf); % time steps 

tolerance_sec=inf;    % Tolerance in sec (if inf --> always include/exclude when overlappping annotations)
tolerance=tolerance_sec*sf;

correction='include'; % or exclude

calculate_agreement=true; % this only makes sense if running the script for multiple subjects

%% loop over subjects
for s=1:length(subjects)
  %% find the annotation files of both raters for this subject
  files{1}=dir(fullfile(folder_rater1, sprintf('**/*%s*.txt', subjects{s})));
  files{2}=dir(fullfile(folder_rater2, sprintf('**/*%s*.txt', subjects{s})));  
  % check if the same number of files are found for the two annotators
  if length(files{1})~= length(files{2})
    error('Not the same number of annotation files (.txt) found in the two annotator folders for subject %s',subjects{s}) 
  else
    % sort names in ascending order so order of files for the two rater
    % corresond (make sure they have similar names though!)
    for i=1:2
      [~, idx]=sort({files{i}.name});
      files{i}=files{i}(idx);
    end
  end
  
  %% loop over files
  for f=1:length(files{1})
    file_rater1=fullfile(files{1}(f).folder, files{1}(f).name);
    file_rater2=fullfile(files{2}(f).folder, files{2}(f).name);
    file_combined=fullfile(folder_combined, sprintf('sub-%s_file-%02d.tsv', subjects{s}, f));
    agreement_table=fullfile(folder_combined, 'agreement_table.tsv');
    combine_FOGannotations(file_rater1, file_rater2, file_combined, agreement_table, subjects{s}, 'include', 2)
  end
end