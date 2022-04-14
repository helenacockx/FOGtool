%% set-up:
clear all; close all;
addpath C:\Users\helen\Documents\PhD\FOG_annotation\matlabscripts\FOG_annotations

% input
folder_rater1='C:\Users\helen\Documents\PhD\FOG_annotation\matlabscripts\FOG_annotations\examplefiles\annotations_rater1';
folder_rater2='C:\Users\helen\Documents\PhD\FOG_annotation\matlabscripts\FOG_annotations\examplefiles\annotations_rater2';
folder_combined='C:\Users\helen\Documents\PhD\FOG_annotation\matlabscripts\FOG_annotations\examplefiles\annotations_combined';
agreement_table=fullfile(folder_combined, 'agreement_table.tsv');
tolerance_sec=2;    % Tolerance in sec
correction='include'; % or exclude

%% find the annotation files of both raters for this subject
files{1}=dir(fullfile(folder_rater1, '**/*.txt'));
files{2}=dir(fullfile(folder_rater2, '**/*.txt'));
% check if the same number of files are found for the two annotators
if length(files{1})~= length(files{2})
  error('Not the same number of annotation files (.txt) found in the two annotator folders')
elseif isempty(files{1}) | isempty(files{2})
  error('No annotation files (.txt) found')
else
  % sort names in ascending order so order of files for the two rater
  % correspond (make sure they have similar names though!)
  for i=1:2
    [~, idx]=sort({files{i}.name});
    files{i}=files{i}(idx);
  end
end
  
%% loop over files
for f=1:length(files{1})
  file_rater1=fullfile(files{1}(f).folder, files{1}(f).name);
  file_rater2=fullfile(files{2}(f).folder, files{2}(f).name);
  for i=1:2
    file_parts{i} = strsplit(files{i}(f).name, '_');
    sub{i} = file_parts{i}{1};
    middle{i} = strjoin(file_parts{i}(2:end-1), '_');
    annotator{i} = file_parts{i}{end};
    full_file{i} = fullfile(files{i}(f).folder, files{i}(f).name);
  end
  if isempty(middle{1})
    name_combined = sub{1};
  else
    name_combined = strcat(sub{1}, '_', middle{1});
  end
  fprintf('<strong>%s</strong>\n', name_combined)
  file_combined=fullfile(folder_combined, name_combined);
  combine_FOGannotations(file_rater1, file_rater2, file_combined, agreement_table, sub{1}, 'include', 2, 1,1);
end

%% calculate final agreement
agr_table = readtable(agreement_table, 'FileType', 'text', 'Delimiter', '\t');
agr_table = agr_table(~strcmp(agr_table.filename, 'TOTAL'),:); % remove the TOTAL row if already present
agreement = agreement_calculator(agr_table); % calculate agreement
m = size(agr_table, 2);

% create new row containing total agreement
last_row = table('Size', [1 m], 'VariableTypes', ['cell' 'cell' repmat({'double'}, 1, 16)], 'VariableNames', agr_table.Properties.VariableNames);
last_row.filename = {'TOTAL'};
last_row.positive_agreement = agreement.pos_agree;
last_row.negative_agreement = agreement.neg_agree;
last_row.prevalence_index= agreement.prev_indx;
w = agr_table.duration_FOG_agreed/sum(agr_table.duration_FOG_agreed); % attribute weigths according to the total duration of agreed FOG
last_row.agreement_trigger = nansum(w.*agr_table.agreement_trigger);
last_row.agreement_type = nansum(w.*agr_table.agreement_type);
last_row{end,8:end} = sum(agr_table{:,8:end});

% show table
final_table = [agr_table; last_row];
fprintf('<strong>FINAL AGREEMENT:</strong>\n')
display(final_table(:,2:end));

% save new table
[path, name, ext] = fileparts(agreement_table);
writetable(final_table, agreement_table,  'FileType', 'text', 'Delimiter', '\t'); % as tsv
writetable(final_table, fullfile(path, [name, '.xlsx']), 'FileType', 'spreadsheet'); % as spread sheet