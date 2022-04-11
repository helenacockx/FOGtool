function [agreement_info] = combine_FOGannotations(filename_rater1, filename_rater2, filename_combined, filename_agreement_table, ID, correction, tolerance_sec, show_agreement, save_fig)
% Script to combine and compare annotations of two raters.

% STEPS:
% 1. export ELAN annotations as Tab-delimited text (seperate column for each tier: true; at least include begin time
% and end time in msec)
% 2. combine annnotations of two raters by running this script. This will
% generate a new .tsv file in the folder_combined that can be imported in ELAN
% 3. import the new .tsv file in ELAN. This will generate a tier
% 'FOG_agreed_Trigger', 'FOG_agreed_Type', 'FOG_disagreed_Trigger',
% 'FOG_disagreed_Type', 'check_annotation', 'NOTES_rater1', 'NOTES_rater2',
% and potential extra tiers, taken over from both raters.
% - create new ELAN file; make sure to add the same video offsets
% - import CSV/tab-delimited text file: the columnc should automatically be 
% named correctly: 'Begin Time' (in msec), 'End Time' (in msec), 'Tier',
% and 'Annotation'; specify first row of data =2; specify delimiter = tab; Skip
% empty cells, don't create empty annotations = true.
% 4. the raters discuss (together with a third rater) which of the annotations
% in FOG_disagreed_Trigger/FOG_disagreed_Type to keep it or delete.
% When an annotation is flagged with 'check_type' or
% 'check_trigger', the raters did not agree on the characterisation of 
% the annoation, and the value of the type or trigger should be checked respectively
% 5. the new file can be saved and exported again to a fileformat of
% choice. By combining the agreed and disagreed FOG's (which have a
% consensus now), you have all the FOG's.

% TO DO:
% - what if 'export multiple files as'

% - check script with empty annotation files
%% Turn on for debugging VS
% tolerance_sec=2;
% filename_agreement_table=agreement_table;
% ID=subjects(VSnmmr2).name;
% correction='include';
%% set-up:
sf=1000;    % choose sampling frequency
ts=(1/sf);  % time steps
tolerance=tolerance_sec*sf;
flag_noGaittask = false;
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames'); % Surpress error notification

flag_nogaittask = false;
flag_notrigger = false;
flag_notype = false;

%% import annotations
opts_rater1 = detectImportOptions(filename_rater1); % read the import options
opts_rater1 = setvartype(opts_rater1,{'FOG_Trigger','FOG_Type'},'char'); % set to char, default is otherwise double
opts_rater2 = detectImportOptions(filename_rater2);        
opts_rater2 = setvartype(opts_rater2,{'FOG_Trigger','FOG_Type'},'char');   
annotations{1}=readtable(filename_rater1, opts_rater1);
annotations{2}=readtable(filename_rater2, opts_rater2);

% convert all variable names to lower case (for ELAN version compatibility)
annotations{1}.Properties.VariableNames = lower(annotations{1}.Properties.VariableNames);
annotations{2}.Properties.VariableNames = lower(annotations{2}.Properties.VariableNames);

% extract FOG annotations and gait tasks
for i=1:2
  
  % REMOVE THESE LINES IN THE DEFINITIVE VERSION!!!  
  % convert time from seconds to msec for annotations Helena
    try
      annotations{i}.begintime_msec = annotations{i}.begintime_ss_msec*1000;
      annotations{i}.endtime_msec = annotations{i}.endtime_ss_msec*1000;
    end
    
    FOG_annotations{i}=annotations{i}(~ismissing(annotations{i}.fog_trigger(:))|~ismissing(annotations{i}.fog_trigger(:)),:);
    FOG_annotations{i} = FOG_annotations{i}(:, {'begintime_msec', 'endtime_msec', 'fog_trigger', 'fog_type'});% remove extra columns
    
    % check if each FOG has been labeled with a FOG_Trigger
    idx = find(ismissing(FOG_annotations{i}.fog_trigger(:)));
    if ~isempty(idx)
      if length(idx) == height(FOG_annotations{i})
        warning('Rater %.0d did not characterize any FOG by trigger', i)
        flag_notrigger = true;
      else 
        warning('Not all FOG events were both annotated for FOG_Trigger for rater %.0d', i)
        display(FOG_annotations{i}(idx,:))
      end
    end
    % check if each FOG has been labeled with a FOG_Type
    idx = find(ismissing(FOG_annotations{i}.fog_type(:)));
    if ~isempty(idx)
      if length(idx) == height(FOG_annotations{i})
        warning('Rater %.0d did not characterize any FOG by type', i)
        flag_notype = true;
      else 
        warning('Not all FOG events were both annotated for FOG_Type for rater %.0d', i)
        display(FOG_annotations{i}(idx,:))
      end
    end
    
    % check if all FOG events are seperated by at least 2 msec
    idx = find(([FOG_annotations{i}.begintime_msec; nan] - [nan; FOG_annotations{i}.endtime_msec])<2);
    FOG_annotations{i}.begintime_msec(idx) = FOG_annotations{i}.begintime_msec(idx)+2;
    
    % convert table from wide to long format
    varnames = annotations{i}.Properties.VariableNames; % find extra tiers names
    annotations_long{i} = stack(annotations{i}, varnames(~contains(varnames, {'begintime', 'endtime', 'duration'})), 'IndexVariableName', 'Tier', 'NewDataVariableName', 'Annotation');
    annotations_long{i} = rmmissing(annotations_long{i});
    annotations_long{i} = annotations_long{i}(:, {'begintime_msec', 'endtime_msec', 'Tier', 'Annotation'});% remove extra columns
    
    % calculate total duration based on the gait tasks
    gait_tasks{i}=annotations_long{i}(annotations_long{i}.Tier == 'gait_task', :);
    if isempty(gait_tasks{i})
        warning('No gait_tasks were found for rater %.0d', i)
        duration_gait_tasks{i}=max(FOG_annotations{i}.endtime_msec)+1000;
    else
        duration_gait_tasks{i}=sum(gait_tasks{i}.endtime_msec-gait_tasks{i}.begintime_msec);
    end
end

% calculate total_duration and endtime for the annotations
if ~isempty(gait_tasks{1}) & ~isempty(gait_tasks{2})
    % check if total_duration is the same for both files
    if round(duration_gait_tasks{1}, -3)~=round(duration_gait_tasks{2}, -3)
        warning('total duration of gait tasks was not the same for both raters. Using the gait_tasks of rater 1.')
        fprintf('total duration of rater 1 (sec): %d \n', round(duration_gait_tasks{1}/1000));
        fprintf('total duration of rater 2 (sec): %d \n', round(duration_gait_tasks{2}/1000));
    end
    total_duration = duration_gait_tasks{1};
    endtime = max(gait_tasks{1}.endtime_msec);
    gait_tasks=gait_tasks{1};
elseif ~isempty(gait_tasks{1}) | ~isempty(gait_tasks{2})
    rater = find([~isempty(gait_tasks{1})  ~isempty(gait_tasks{2})]);
    warning('Only annotations of rater %.0d contained gait_tasks. Using those to calculate total duration.', rater)
    total_duration = duration_gait_tasks{rater};
    endtime = max(gait_tasks{rater}.endtime_msec);
    gait_tasks=gait_tasks{rater};
elseif isempty(gait_tasks{1}) & isempty(gait_tasks{2})
    warning('No gait tasks were found for any of the raters. Assuming that the gait_task ended 1 sec after the last FOG event. The agreement parameters will not be calculated and the visualization might be wrong.')
    flag_nogaittask = true;
    total_duration = max([duration_gait_tasks{1}, duration_gait_tasks{2}]);
    endtime = total_duration;
    if isempty(endtime)
        endtime = 0;
        total_duration = 0;
    end
    gait_tasks=table(0, endtime, 'VariableNames', {'begintime_msec', 'endtime_msec'});
end

% convert annotations to boolean vectors based on the given sampling frequency
t=(0:1:endtime+1)*ts; % time vector (add 1 extra seconds to make sure that the boolvec_FOG goes back to zero after the last FOG)
boolvec_task=nan(1,(length(t))); % boolean vector including all time points
% make boolvec 0 during gait_tasks
for i=1:height(gait_tasks)
    boolvec_task(round(gait_tasks.begintime_msec(i))+1:round(gait_tasks.endtime_msec(i))+1)=0; % +1 because going from msec to samples
end

for i=1:2
    % add extra column with begin sample and end sample of FOG annotation, nu
    % overbodig?
    % FOG_annotations{i}.begintime_msec=round(FOG_annotations{i}.begintime_Ss_msec*sf+1);
    % FOG_annotations{i}.endtime_msec=round(FOG_annotations{i}.endtime_Ss_msec*sf+1);
    
    % create a boolean vector with the FOG annotations of this rater
    boolvec_FOG=boolvec_task;
    for k=1:height(FOG_annotations{i})
        boolvec_FOG((round(FOG_annotations{i}.begintime_msec(k))+1):(round(FOG_annotations{i}.endtime_msec(k))+1))=1; % +1 because going from msec to samples
    end
    % check if all FOGs are falling inside the gait_task
    if sum(boolvec_FOG ==1 & isnan(boolvec_task))>0
      startsample = find(diff([boolvec_FOG ==1 & isnan(boolvec_task)])==1);
      idx = [];
      for k=1:length(startsample)
        idx = [idx; find(startsample(k)>=FOG_annotations{i}.begintime_msec & startsample(k)<=FOG_annotations{i}.endtime_msec)];
      end
      warning('%d FOG events are falling outside the gait_task and are trimmed or removed from the list:', numel(startsample))
      display(FOG_annotations{i}(idx,:))
      display(gait_tasks)
    end
    FOG_vector{i}=boolvec_FOG + boolvec_task; % + boolvec_task to make sure that FOGs falling outside the gait_task are made nan
    % make sure all FOG events go back to 0 before gait_task ends and
    % starts from 0 when gait_task starts
    probl_edge = find(FOG_vector{i}==1 & (isnan([diff(FOG_vector{i}) nan])|isnan([nan diff(FOG_vector{i})]))); % find problematic edges (= [1 nan] or [nan 1])
    FOG_vector{i}(probl_edge) = 0;
end

%% combine annotations
% 0 = definitely no FOG (white); 2 = definitely FOG (black); 1 = possible FOG (non-overlapping, grey area); nan = no gait_task
FOG_summed=FOG_vector{1}+FOG_vector{2}; 

% find begin and end samples of the possible FOG events
beginsample=find(FOG_summed==1 & diff([0 FOG_summed])~=0);% find beginsamples when going from definetly (black/white) to possible (grey)
endsample=find(FOG_summed==1 &  diff([FOG_summed 0])~=0); % find endsamples when going back from possible (grey) to definitely (black/white)
if length(beginsample)~= length(endsample)
    error('Not all begin and end samples were found for the possible FOG events')
end

% include or exclude possible FOG events based on the chosen parameters
% (correction & tolerance)
FOG_corrected=FOG_summed;
for k=1:length(beginsample)
  % if isolated possible FOG episode(only one rater annotated the episode)
  if FOG_summed(beginsample(k)-1)== 0 & FOG_summed(endsample(k)+1)==0
    FOG_corrected(beginsample(k):endsample(k))=1; % remains grey for discussion
  % if non-isolated possible FOG with a duration > tolerance
  elseif endsample(k)-beginsample(k)> tolerance
    FOG_corrected(beginsample(k):endsample(k))=1; % remains grey for discussion
  % if non-isolated possible FOG with duration < tolerance
  else
    % outcome depends on the correction parameter
    switch correction
      case 'include'
        FOG_corrected(beginsample(k):endsample(k)) = 2; % definitely FOG
      case 'exclude'
        FOG_corrected(beginsample(k):endsample(k)) = 0; % definitely no FOG
      otherwise
        error('no valid correction option was chosen. Please choose include or exclude.')
    end
  end
end

% find the agreed (definitely FOG/no FOG) and disagreed (to-be-discussed) FOG's
FOG_agreed = (FOG_corrected==2) + boolvec_task;
FOG_disagreed = (FOG_corrected==1) + boolvec_task;

%% convert FOG disagreed to a table and add the rater, trigger and type for this FOG
[beginsample, endsample]=vec2event(FOG_disagreed);% find beginsample and endsample of each event
n=length(beginsample);

% pre-allocate output
clear varnames vartypes
varnames = {'begintime_msec', 'endtime_msec','Tier', 'Annotation', 'rater'};
vartypes(1,[1:2 5])={'double'};
vartypes(1,[3:4])={'string'};
FOG_disagreed_t=table('Size', [2*n, length(vartypes)], 'VariableNames', varnames, 'VariableTypes', vartypes); % 2*n because one row FOG_Trigger, and one row FOG_type
FOG_disagreed_t.Tier = repmat({'FOG_disagreed_Trigger'; 'FOG_disagreed_Type'}, n, 1);
FOG_disagreed_t.begintime_msec=repelem(beginsample'-1, 2,1); % 2 rows for each FOG episode
FOG_disagreed_t.endtime_msec=repelem(endsample'-1,2,1);

% find rater, trigger and type of this FOG
for k=1:n
    % find annotations that fall within this event
    idx_rater1 = overlappingevt(FOG_annotations{1}, beginsample(k), endsample(k));
    idx_rater2 = overlappingevt(FOG_annotations{2}, beginsample(k), endsample(k));
    
    if length(idx_rater1)==1 & isempty(idx_rater2) % this is an annotation from rater 1       
        FOG_disagreed_t.rater(2*k-1:2*k)=1;
        FOG_disagreed_t.Annotation(2*k-1)=FOG_annotations{1}.fog_trigger(idx_rater1);
        FOG_disagreed_t.Annotation(2*k)=FOG_annotations{1}.fog_type(idx_rater1);
    elseif isempty(idx_rater1) & length(idx_rater2)==1 % this is an annotation from rater 2
        FOG_disagreed_t.rater(2*k-1:2*k)=2;
        FOG_disagreed_t.Annotation(2*k-1)=FOG_annotations{2}.fog_trigger(idx_rater2);
        FOG_disagreed_t.Annotation(2*k)=FOG_annotations{2}.fog_type(idx_rater2);
    else
        error('Multiple annotations were found for this disagreed FOG episode') 
    end
end

%% convert FOG agreed to a table and check whether trigger and type for this FOG of both raters was the same (if not combine both values and check_trigger/check_type=true)
[beginsample, endsample]=vec2event(FOG_agreed);% find beginsample and endsample of each event
n=length(beginsample);
FOG_agreed_t=table('Size', [2*n, length(vartypes)], 'VariableNames', varnames, 'VariableTypes', vartypes); % 2*n because one row FOG_Trigger, and one row FOG_type
FOG_agreed_t.Tier = repmat({'FOG_agreed_Trigger'; 'FOG_agreed_Type'}, n, 1);
FOG_agreed_t.begintime_msec=repelem(beginsample'-1, 2,1); % 2 rows for each FOG episode
FOG_agreed_t.endtime_msec=repelem(endsample'-1,2,1);

for k=1:n
    % find annotations of the two raters that fall within this event
    idx_rater1 = overlappingevt(FOG_annotations{1}, beginsample(k), endsample(k));
    idx_rater2 = overlappingevt(FOG_annotations{2}, beginsample(k), endsample(k));
    if isempty(idx_rater1) | isempty(idx_rater2)
        error('No annotations were found for this agreed FOG episode')
    end

    % check trigger
    triggers=[FOG_annotations{1}.fog_trigger(idx_rater1); FOG_annotations{2}.fog_trigger(idx_rater2)];
    triggers=unique(triggers); % unique values for triggers
    if length(triggers)==1 | flag_notrigger % the same value was given for FOG_Trigger/only one annotator characterized the trigger 
        FOG_agreed_t.Annotation(2*k-1)=triggers(~ismissing(triggers));
    else % a different value was given for trigger
        % make extra annotation to check trigger
        FOG_agreed_t(end+1,:) = FOG_agreed_t(2*k-1,:);
        FOG_agreed_t.Tier(end) = {'check_annotation'};
        FOG_agreed_t.Annotation(end)={'check_trigger'};
        % combine both values to one string
        triggers(cellfun(@isempty, triggers)) = {'not specified'}; % replace empty cell, by 'not specified'
        trig_combi=triggers{1};
        for j=2:length(triggers)
            trig_combi=[trig_combi ' / ' triggers{j}];
        end
        FOG_agreed_t.Annotation(2*k-1)={trig_combi};
    end

    % check type
    types=[FOG_annotations{1}.fog_type(idx_rater1); FOG_annotations{2}.fog_type(idx_rater2)];
    types=unique(types);% unique values for types
    if length(types)==1 | flag_notype % the same value was given for FOG_Type/only one annotator characterized the type
        FOG_agreed_t.Annotation(2*k)=types(~ismissing(types));
    else  % a different value was given for type
        % make extra annotation to check type
        FOG_agreed_t(end+1,:) = FOG_agreed_t(2*k,:);
        FOG_agreed_t.Tier(end) = {'check_annotation'};
        FOG_agreed_t.Annotation(end)={'check_type'};
        % combine both values to one string
        types(cellfun(@isempty, types)) = {'not specified'}; % replace empty cell, by 'not specified'
        type_combi=types{1};
        for h=2:length(types)
            type_combi=[type_combi ' / ' types{h}];
        end
        FOG_agreed_t.Annotation(2*k)={type_combi};
    end
end


%% Visualize the results 
image.save = save_fig;
image.name = filename_combined;
PlotAnn(FOG_vector, FOG_agreed, FOG_disagreed, gait_tasks, t, image)

%% combine the agreed and disagreed tables into one table and extra tiers
% combine agreed and disagreed FOG episodes
FOG_all_t=[FOG_agreed_t; FOG_disagreed_t]; 
FOG_all_t=removevars(FOG_all_t, 'rater'); % remove rater from table, so discussion in ELAN is blind for this

 % add gait_tasks
final_table = [FOG_all_t; gait_tasks]; 

% add notes and extra tiers 
all_extra_tiers = [];
for i=1:2 
  notes = annotations_long{i}(annotations_long{i}.Tier=='notes',:);
  if i==1
    notes.Tier = categorical(repmat("NOTES_rater1", height(notes),1));
  elseif i==2
    notes.Tier = categorical(repmat("NOTES_rater2", height(notes),1));
  end
  extra_tiers = annotations_long{i}(all(annotations_long{i}.Tier ~= {'fog_trigger', 'fog_type', 'gait_task', 'notes'},2),:);
  all_extra_tiers = [all_extra_tiers; notes; extra_tiers];
end
% remove duplicate rows
all_extra_tiers = sortrows(all_extra_tiers);
all_extra_tiers_rounded = table(round(all_extra_tiers.begintime_msec, -2), round(all_extra_tiers.endtime_msec, -2), all_extra_tiers.Annotation); % annotations can differ by max 100 msec
[~, unique_idx] = unique(all_extra_tiers_rounded, 'rows');
all_extra_tiers = all_extra_tiers(unique_idx,:);
final_table = [final_table; all_extra_tiers];

% save table
header_names = {'Begin Time', 'End Time', 'Tier', 'Annotation'};
final_table_cell = [header_names; table2cell(final_table)];
writecell(final_table_cell, strcat(filename_combined, '_annotations-combined'), 'Filetype', 'text', 'Delimiter', '\t'); % workaround to add spaces in header names

%% fill in agreement table
% create new vector that also contains information about the rater of the
% disagreed FOG
FOG_summed_v2=FOG_vector{1}+2*FOG_vector{2}; % 0=agreed no FOG; 3=agreed FOG; 1=FOG only annotated by rater 1; 2=FOG only annotated by rater 2

% make an agreement table for this file
varnames={'subject', 'filename', 'positive_agreement', 'negative_agreement', 'prevalence_index',...
    'agreement_trigger', 'agreement_type', 'number_FOG_rater1', 'number_FOG_rater2', 'number_FOG_agreed',...
    'number_FOG_disagreed_rater1', 'number_FOG_disagreed_rater2',...
    'duration_FOG_rater1', 'duration_FOG_rater2', 'duration_FOG_agreed', ...
    'duration_FOG_disagreed_rater1', 'duration_FOG_disagreed_rater2', 'total_duration'}; % 'kappa', 'ICC',
vartypes=[repmat({'string'}, [1,2]), repmat({'double'}, [1,16])];
agreement_info=table('Size', [1, 18], 'VariableNames', varnames, 'VariableTypes', vartypes);


agreement_info.subject={ID};
[path, name, ext]=fileparts(filename_combined);
agreement_info.filename={name};
agreement_info.number_FOG_rater1=height(FOG_annotations{1});
agreement_info.duration_FOG_rater1=sum([FOG_annotations{1}.endtime_msec-FOG_annotations{1}.begintime_msec])/1000;
agreement_info.number_FOG_rater2=height(FOG_annotations{2});
agreement_info.duration_FOG_rater2=sum([FOG_annotations{2}.endtime_msec-FOG_annotations{2}.begintime_msec])/1000;
agreement_info.number_FOG_agreed=sum(strcmp(FOG_agreed_t.Tier, 'FOG_agreed_Trigger')); % only for info, not to calculate agreement (because uses adjusted FOG annotations)
agreement_info.duration_FOG_agreed=sum(FOG_summed_v2==3)/1000; % in sec
agreement_info.number_FOG_disagreed_rater1=sum(FOG_disagreed_t.rater==1)/2; % only for info, not to calculate agreement (because uses adjusted FOG annotations)
agreement_info.duration_FOG_disagreed_rater1=sum(FOG_summed_v2==1)/1000;
agreement_info.number_FOG_disagreed_rater2=sum(FOG_disagreed_t.rater==2)/2; % only for info, not to calculate agreement (because uses adjusted FOG annotations)
agreement_info.duration_FOG_disagreed_rater2=sum(FOG_summed_v2==2)/1000;
agreement_info.total_duration=total_duration/1000;

[agreement] = agreement_calculator(agreement_info);
agreement_info.positive_agreement = agreement.pos_agree;
agreement_info.negative_agreement = agreement.neg_agree;
agreement_info.prevalence_index = agreement.prev_indx;
  
if flag_nogaittask
  warning('No gait tasks were found, do not calculate agreement parameters')
  agreement_info.total_duration = nan;
  agreement_info.positive_agreement = nan;
  agreement_info.negative_agreement = nan;
  agreement_info.prevalence_index = nan;
end

% calculate %agreement on trigger and type
if flag_notrigger
  agreement_info.agreement_trigger = nan;
else
  agreement_info.agreement_trigger = (agreement_info.number_FOG_agreed-sum(strcmp(FOG_agreed_t.Annotation, 'check_trigger')))/agreement_info.number_FOG_agreed;
end
if flag_notype
  agreement_info.agreement_type = nan;
else 
  agreement_info.agreement_type = (agreement_info.number_FOG_agreed-sum(strcmp(FOG_agreed_t.Annotation, 'check_type')))/agreement_info.number_FOG_agreed;
end

% display
if show_agreement
  fprintf('Agreement info of file %s: \n', name)
  display(table(varnames(3:end)', round(agreement_info{:,3:end}',2), 'VariableNames', {'annotation_info', 'value'}))
end

% load the big agreement table if present
if exist(filename_agreement_table, 'file')
    agreement_t=readtable(filename_agreement_table, 'FileType', 'text', 'ReadVariableNames', 1, 'HeaderLines', 0);
    % check if this file is already part of the agreement table
    if any(strcmp(agreement_t.filename, name)) % update the info in the table
        n=find(strcmp(agreement_t.filename, name));
        agreement_t(n,:)=agreement_info;
    else % add an extra row to the table
        agreement_t=[agreement_t; agreement_info];
    end
else % create a new table
    agreement_t=agreement_info;
end

% save the table
writetable(agreement_t, filename_agreement_table,  'FileType', 'text', 'Delimiter', '\t');


%% HELPER FUNCTIONS
% VEC2EVENT
function     [beginsample, endsample]=vec2event(boolvec)
tmp = diff([0 boolvec 0]);
beginsample = find(tmp==+1);
endsample = find(tmp==-1) - 1;

% PlotAnn
function PlotAnn(FOG_vector, FOG_agreed, FOG_disagreed, gait_tasks, t, image)

% open new figure
ImageResult = figure();
set(gca, 'color', [0.9 0.9 0.9])
color_discuss = sscanf('9b9b9b','%2x%2x%2x',[1 3])/255;

% plot gait_tasks
hold on
for o=1:height(gait_tasks)
  rectangle('Position',[t(gait_tasks.begintime_msec(o)+1) 0 ...
    t(gait_tasks.endtime_msec(o)-gait_tasks.begintime_msec(o)) 1],...
    'linestyle', '-', 'FaceColor', 'w')
  rectangle('Position',[t(gait_tasks.begintime_msec(o)+1) 2 ...
    t(gait_tasks.endtime_msec(o)-gait_tasks.begintime_msec(o)) 1],...
    'linestyle', '-', 'FaceColor', 'w')
  rectangle('Position',[t(gait_tasks.begintime_msec(o)+1) 4 ...
    t(gait_tasks.endtime_msec(o)-gait_tasks.begintime_msec(o)) 1],...
    'linestyle', '-', 'FaceColor', 'w')
end

% plot FOG rater 1
[beginsample_FOG_rater1, endsample_FOG_rater1]=vec2event(FOG_vector{1});
for o=1:length(beginsample_FOG_rater1)
  rectangle('Position', [t(beginsample_FOG_rater1(o)) 4 ...
    t(endsample_FOG_rater1(o)-beginsample_FOG_rater1(o)) 1],...
    'FaceColor', 'k')
end

% plot FOG rater 2
[beginsample_FOG_rater2, endsample_FOG_rater2]=vec2event(FOG_vector{2});
for o=1:length(beginsample_FOG_rater2)
  rectangle('Position', [t(beginsample_FOG_rater2(o)) 2 ...
    t(endsample_FOG_rater2(o)-beginsample_FOG_rater2(o)) 1],...
    'FaceColor', 'k')
end

% plot FOG agreed
[beginsample_FOGagreed, endsample_FOGagreed]=vec2event(FOG_agreed);
for o=1:length(beginsample_FOGagreed)
  rectangle('Position', [t(beginsample_FOGagreed(o)) 0 ...
    t(endsample_FOGagreed(o)-beginsample_FOGagreed(o)) 1],...
    'FaceColor', 'k')
end

% plot remaining grey areas
[beginsample_FOGdisagreed, endsample_FOGdisagreed]=vec2event(FOG_disagreed);
for o=1:length(beginsample_FOGdisagreed)
  rectangle('Position', [t(beginsample_FOGdisagreed(o)) 0 ...
    t(endsample_FOGdisagreed(o)-beginsample_FOGdisagreed(o)) 1],...
    'FaceColor', color_discuss, 'EdgeColor', 'none')
end

% handle axes and add labels 
ylim([-1 6])
yticks([0.5 2.5 4.5])
yticklabels({'outcome', 'rater 2', 'rater 1'})
xlabel('time (in s)')

% add legend 
dummy_FOG = fill(0,0, 'k');
dummy_noFOG = fill(0,0, 'w');
dummy_discuss = fill(0, 0, color_discuss);
legend([dummy_FOG, dummy_noFOG, dummy_discuss], {'FOG', 'no FOG', 'discuss'}, ...
  'orientation', 'horizontal', 'Location', 'north')

% format layout
set( findall(ImageResult, '-property', 'fontsize'), 'fontsize', 11);
% ImageResult.WindowState = 'maximized';
[path, name, ext]=fileparts(image.name);
title(strrep(name, '_', ' '),'fontsize', 18);

% save image
if image.save == 1
  saveas(ImageResult,[image.name '.png']);
end

% OVERLAPPINGEVT
function [idx] = overlappingevt(annotations, beginsample, endsample)
% find the indices of annotation events that fall within the event with the
% given [beginsample endsample].
idx=find(([annotations.begintime_msec]<=beginsample & [annotations.endtime_msec]>beginsample) |... % annotation includes the beginsample
    ([annotations.begintime_msec]<endsample & [annotations.endtime_msec]>=endsample) | ... % annotation includes the endsample
    ([annotations.begintime_msec]>=beginsample & [annotations.endtime_msec]<endsample)); % annotation falls within the event
