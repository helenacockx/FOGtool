function combine_FOGannotations(filename_rater1, filename_rater2, filename_combined, filename_agreement_table, ID, correction, tolerance_sec)
% Script to combine and compare annotations of two raters.

% STEPS:
% 1. export ELAN annotations as Tab-delimited text (seperate column for each tier: true; at least include begin time
% and end time in ss.msec)
% 2. combine annnotations of two raters by running this script. This will
% generate a new .tsv file in the folder_combined that can be imported in ELAN
% 3. import the new .tsv file in ELAN. This will generate a tier
% 'FOG_agreed_Trigger', 'FOG_agreed_Type', 'FOG_disagreed_Trigger',
% 'FOG_disagreed_Type', 'check_trigger', and 'check_type.
% 4. the third rater decides on all the annotations in FOG_disagreed_... to
% keep it or delete it. When an annotation is flagged with 'check_type' or
% 'check_trigger', the raters did not agree on the value of the annoation,
% and the value of the type or trigger should be checked respectively
% 5. the new file can be saved and exported again to a fileformat of
% choice. By combining the agreed and disagreed FOG's (which have a
% consensus now), you have all the FOG's. 

% TO DO:
% - what if 'export multiple files as'
% - what if FOG_Trigger and FOG_Type are not dependent tiers?/are not both
% in the same line
% - option to visualize/always visualize? use time vector that corresponds
% to ELAN? but also depends on offset master media file (so uncheck)
% - check script with empty annotation files

%% set-up:
sf=1000; % choose sampling frequency
ts=(1/sf); % time steps 
tolerance=tolerance_sec*sf;  
 warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

%% collect the annotations of these files
annotations{1}=readtable(filename_rater1, 'ReadVariableNames', 1, 'HeaderLines', 0);
annotations{2}=readtable(filename_rater2, 'ReadVariableNames', 1, 'HeaderLines', 0);
for i=1:2
  % only select FOG annotations & notes
  FOG_annotations{i}=annotations{i}(~ismissing(annotations{i}.FOG_Trigger(:))|~ismissing(annotations{i}.FOG_Trigger(:)),:);
  NOTES{i}=annotations{i}(~ismissing(annotations{i}.NOTES(:)) & (ismissing(annotations{i}.FOG_Trigger(:))&ismissing(annotations{i}.FOG_Type(:))),:); % extra notes that were not annotated during a FOG
  % check if each FOG has been labeled with a FOG_Trigger and a FOG_Type
  idx=find(ismissing(FOG_annotations{i}.FOG_Trigger(:))|ismissing(FOG_annotations{i}.FOG_Type));
  if ~isempty(idx)
    warning('Not all FOG events were both annotated for FOG_Trigger and FOG_Type for rater %.0d', i)
    display(FOG_annotations{i}(idx,:))
  end
  % calculate total duration based on the gait tasks
  gait_tasks{i}=annotations{i}(~ismissing(annotations{i}.gait_task(:)),:);
  if isempty(gait_tasks{i})
    warning('No gait_tasks were found for rater %.0d. Assuming that the gait_task ended after the last FOG event.', i)
    duration_gait_tasks{i}=max(FOG_annotations{i}.EndTime_Ss_msec);
  else
    duration_gait_tasks{i}=sum(gait_tasks{i}.EndTime_Ss_msec-gait_tasks{i}.BeginTime_Ss_msec);
  end
end

% calculate total_duration and endtime for the annotations
if ~isempty(gait_tasks{1}) & ~isempty(gait_tasks{2})
  % check if total_duration is the same for both files
  if round(duration_gait_tasks{1})~=round(duration_gait_tasks{2})
    warning('total duration of gait tasks was not the same for both raters. Using the gait_tasks of rater 1.')
    fprintf('total duration of rater 1: %d \n', round(duration_gait_tasks{1}));
    fprintf('total duration of rater 2: %d \n', round(duration_gait_tasks{2}));
  end
  total_duration = duration_gait_tasks{1};
  endtime = max(gait_tasks{1}.EndTime_Ss_msec);
  gait_tasks=gait_tasks{1};
elseif ~isempty(gait_tasks{1}) | ~isempty(gait_tasks{2})
  rater = find([~isempty(gait_tasks{1})  ~isempty(gait_tasks{2})]);
  warning('Only annotations of rater %.0d contained gait_tasks. Using those to calculate total duration', rater)
  total_duration = duration_gait_tasks{rater};
  endtime = max(gait_tasks{rater}.EndTime_Ss_msec);
  gait_tasks=gait_tasks{rater};
elseif isempty(gait_tasks{1}) & isempty(gait_tasks{2})
  warning('No gait tasks were found for any of the raters. Assuming that the gait_task ended after the last FOG event. Be carefull with the interpretation of the agreement coefficients.')
  total_duration = max([duration_gait_tasks{1}, duration_gait_tasks{2}]);
  endtime = total_duration;
  if isempty(endtime)
    endtime = 0;
    total_duration = 0;
  end
  gait_tasks=table(0, endtime, 'VariableNames', {'BeginTime_Ss_msec', 'EndTime_Ss_msec'});
end

%% convert annotations to boolean vectors based on the given sampling frequency
t=[0:1/sf:endtime+1]; % time vector (add 1 extra seconds to make sure that the boolvec_FOG goes back to zero after the last FOG)
boolvec_task=nan(1,length(t)); % boolean vector including all time points
% add extra column with begin sample and end sample of gait tasks
gait_tasks.BeginTime_sample = round(gait_tasks.BeginTime_Ss_msec*sf + 1);
gait_tasks.EndTime_sample = round(gait_tasks.EndTime_Ss_msec*sf + 1);
% make boolvec 0 during gait_tasks
for i=1:height(gait_tasks)
  boolvec_task(gait_tasks.BeginTime_sample(i):gait_tasks.EndTime_sample(i)) = 0; 
end
for i=1:2
  % add extra column with begin sample and end sample of FOG annotation
  FOG_annotations{i}.BeginTime_sample=round(FOG_annotations{i}.BeginTime_Ss_msec*sf+1);
  FOG_annotations{i}.EndTime_sample=round(FOG_annotations{i}.EndTime_Ss_msec*sf+1);
  % create a boolean vector with the FOG annotations of this rater
  boolvec_FOG=boolvec_task;
  for k=1:height(FOG_annotations{i})
    boolvec_FOG(FOG_annotations{i}.BeginTime_sample(k):FOG_annotations{i}.EndTime_sample(k))=1;
  end
  FOG_vector{i}=boolvec_FOG + boolvec_task; % + boolvec_task to make sure that FOGs falling outside the gait_task are made nan
  % make sure all FOG events go back to 0 before gait_task ends
  probl_edge = find(FOG_vector{i}==1 & (isnan([diff(FOG_vector{i}) nan])|isnan([nan diff(FOG_vector{i})]))); % find problematic edges (= [1 nan] or [nan 1])
  FOG_vector{i}(probl_edge) = 0;
end

%% combine annotations
FOG_summed=FOG_vector{1}+FOG_vector{2}; % 0=agreed no FOG; 2=agreed FOG; 1=disagreed FOG; nan=no gait_task

% find begin and end samples of disagreed events
beginsample=find(FOG_summed==1 & diff([0 FOG_summed])~=0);% find beginsamples when going from agreed to disagreed
endsample=find(FOG_summed==1 &  diff([FOG_summed 0])~=0); % find endsamples when going back from disagreed to agreed
if length(beginsample)~= length(endsample)
  error
end

% include or exclude disagreed events based on the chosen parameters
% (correction & tolerance)
FOG_corrected=FOG_summed;
switch correction
  case 'include'
    x=2;
  case 'exclude'
    x=0;
  otherwise
    error('no valid correction option was chosen. Please choose include or exclude.')
end
for k=1:length(beginsample)
  % disagreement on annotation falls within another annotation:
  if FOG_summed(beginsample(k)-1)== 2 & FOG_summed(endsample(k)+1)==2
    if beginsample(k)-endsample(k)<= tolerance
      FOG_corrected(beginsample(k):endsample(k))=x; % or always include?
    else
      FOG_corrected(beginsample(k):endsample(k))=1;
    end
  % only one annotater annotated this FOG episode:
  elseif FOG_summed(beginsample(k)-1)== 0 & FOG_summed(endsample(k)+1)==0
    if endsample(k)-beginsample(k)<= tolerance
      FOG_corrected(beginsample(k):endsample(k))=1; % not agreed annotation 
    else
      FOG_corrected(beginsample(k):endsample(k))=1; % not agreed annotation
    end
  % one annotator annotated the FOG episode earlier than the other
  % annotator:
  elseif FOG_corrected(beginsample(k)-1)== 0 & FOG_summed(endsample(k)+1)==2
    if endsample(k)-beginsample(k)<= tolerance %0.5*
      FOG_corrected(beginsample(k):endsample(k))=x;
    else
      FOG_corrected(beginsample(k):endsample(k))=1; % not agreed annotation
    end
  % one annotator annotated the FOG epidosde longer than the other
  % annotator
  elseif FOG_corrected(beginsample(k)-1)== 2 & FOG_summed(endsample(k)+1)==0
    if endsample(k)-beginsample(k)<= tolerance %0.5*
      FOG_corrected(beginsample(k):endsample(k))=x;
    else
      FOG_corrected(beginsample(k):endsample(k))=1; % not agreed annotation
    end
  end
end

% find the agreed and disagreed FOG's
FOG_agreed = (FOG_corrected==2) + boolvec_task;
FOG_disagreed = (FOG_corrected==1) + boolvec_task;
    
%% visualize
[path, name, ext]=fileparts(filename_combined);
figure;
ax(1)=subplot(4,1,1); plot(t, FOG_vector{1}, 'Color', '#EDB120'); ylim([-1 2]); ylabel('FOG rater 1');title(name);
ax(2)=subplot(4,1,2); plot(t, FOG_vector{2},  'Color','#0072BD'); ylim([-1 2]);ylabel('FOG rater 2'); 
ax(3)=subplot(4,1,3); plot(t, FOG_agreed,  'Color','#77AC30'); ylim([-1 2]);ylabel('FOG agreed')
ax(4)=subplot(4,1,4); plot(t, FOG_disagreed, 'Color', '#A2142F'); ylim([-1 2]);ylabel('FOG disagreed')
linkaxes(ax)
xlabel('time (in sec)');

    
%% convert FOG disagreed to a table and add the rater, trigger and type for this FOG
[beginsample, endsample]=vec2event(FOG_disagreed);% find beginsample and endsample of each event
n=length(beginsample);
FOG_disagreed_t=table(beginsample', endsample',nan(n,1), cell(n,1), cell(n,1),cell(n,1), cell(n,1), cell(n,1), cell(n,1),cell(n,1), cell(n,1), 'VariableNames', {'BeginTime_sample', 'EndTime_sample','rater', 'FOG_agreed_Trigger', 'FOG_agreed_Type', 'FOG_disagreed_Trigger', 'FOG_disagreed_Type', 'check_trigger', 'check_type','NOTES_rater1', 'NOTES_rater2'});% convert into table
for k=1:height(FOG_disagreed_t)
  % find annotations that fall within this event
  idx_rater1 = overlappingevt(FOG_annotations{1}, beginsample(k), endsample(k));
  idx_rater2 = overlappingevt(FOG_annotations{2}, beginsample(k), endsample(k));
  if length(idx_rater1)==1 & isempty(idx_rater2)
    FOG_disagreed_t.rater(k)=1;
    FOG_disagreed_t.FOG_disagreed_Trigger(k)=FOG_annotations{1}.FOG_Trigger(idx_rater1);
    FOG_disagreed_t.FOG_disagreed_Type(k)=FOG_annotations{1}.FOG_Type(idx_rater1);
    try
      FOG_disagreed_t.NOTES_rater1(k)=FOG_annotations{1}.NOTES(idx_rater1);
    end
  elseif isempty(idx_rater1) & length(idx_rater2)==1
    FOG_disagreed_t.rater(k)=2;
    FOG_disagreed_t.FOG_disagreed_Trigger(k)=FOG_annotations{2}.FOG_Trigger(idx_rater2);
    FOG_disagreed_t.FOG_disagreed_Type(k)=FOG_annotations{2}.FOG_Type(idx_rater2);
    try
      FOG_disagreed_t.NOTES_rater2(k)=FOG_annotations{2}.NOTES(idx_rater2);
    end
  else
    error % the disagreed annotations should only be rated by one person and no the event should only overlap with one other event
  end
end

%% convert FOG agreed to a table and check whether trigger and type for this FOG of both raters was the same (if not combine both values and check_trigger/check_type=true)
[beginsample, endsample]=vec2event(FOG_agreed);% find beginsample and endsample of each event
n=length(beginsample);
FOG_agreed_t=table(beginsample', endsample', nan(n,1), cell(n,1), cell(n,1), cell(n,1), cell(n,1), cell(n,1), cell(n,1),cell(n,1), cell(n,1), 'VariableNames', {'BeginTime_sample', 'EndTime_sample', 'rater', 'FOG_agreed_Trigger', 'FOG_agreed_Type', 'FOG_disagreed_Trigger', 'FOG_disagreed_Type',  'check_trigger', 'check_type', 'NOTES_rater1', 'NOTES_rater2'});% convert into table
for k=1:height(FOG_agreed_t)
  % find annotations of the two raters that fall within this event
  idx_rater1 = overlappingevt(FOG_annotations{1}, beginsample(k), endsample(k));
  idx_rater2 = overlappingevt(FOG_annotations{2}, beginsample(k), endsample(k));
  if isempty(idx_rater1) | isempty(idx_rater2)
    error
  end
  triggers=[FOG_annotations{1}.FOG_Trigger(idx_rater1); FOG_annotations{2}.FOG_Trigger(idx_rater2)];
  triggers=unique(triggers); % unique values for triggers
  if length(triggers)==1 % the same value was given for FOG_Trigger
    FOG_agreed_t.FOG_agreed_Trigger(k)=triggers;
  else % a different value was given for trigger
    FOG_agreed_t.check_trigger(k)={'check_trigger'};
    % combine both values to one string
    trig_combi=triggers{1};
    for t=2:length(triggers)
      trig_combi=[trig_combi ' / ' triggers{t}];
    end
    FOG_agreed_t.FOG_agreed_Trigger(k)={trig_combi};
  end
  types=[FOG_annotations{1}.FOG_Type(idx_rater1); FOG_annotations{2}.FOG_Type(idx_rater2)];
  types=unique(types);% unique values for types
  if length(types)==1 % the same value was given for FOG_Type
    FOG_agreed_t.FOG_agreed_Type(k)=types;
  else  % a different value was given for type
    FOG_agreed_t.check_type(k)={'check_type'};
    % combine both values to one string
    type_combi=types{1};
    for t=2:length(types)
      type_combi=[type_combi ' / ' types{t}];
    end
    FOG_agreed_t.FOG_agreed_Type(k)={type_combi};
  end
  try
    FOG_agreed_t.NOTES_rater1(k)=FOG_annotations{1}.NOTES(idx_rater1);
  end
  try
    FOG_agreed_t.NOTES_rater2(k)=FOG_annotations{2}.NOTES(idx_rater2);
  end
end

%% combine the agreed and disagreed tables into one table
% add timing in seconds (instead of samples)
FOG_agreed_t.BeginTime_Ss_msec=(FOG_agreed_t.BeginTime_sample(:)-1)/sf;
FOG_agreed_t.EndTime_Ss_msec=(FOG_agreed_t.EndTime_sample(:)-1)/sf;
FOG_disagreed_t.BeginTime_Ss_msec=(FOG_disagreed_t.BeginTime_sample(:)-1)/sf;
FOG_disagreed_t.EndTime_Ss_msec=(FOG_disagreed_t.EndTime_sample(:)-1)/sf;
FOG_all_t=[FOG_agreed_t; FOG_disagreed_t];

% add extra notes
for k=1:height(NOTES{1})
  FOG_all_t(end+1,:)=[repmat({NaN}, 1,3), repmat({[]}, 1,6) NOTES{1}.NOTES(k) {[]} NOTES{1}.BeginTime_Ss_msec(k) NOTES{1}.EndTime_Ss_msec(k)];
end
for k=1:height(NOTES{2})
  FOG_all_t(end+1,:)=[repmat({NaN}, 1,3), repmat({[]}, 1,7) NOTES{2}.NOTES(k) NOTES{2}.BeginTime_Ss_msec(k) NOTES{2}.EndTime_Ss_msec(k)];
end

% export the table
writetable(FOG_all_t, filename_combined,  'FileType', 'text', 'Delimiter', '\t')


%% fill in agreement table
% create new vector that also contains information about the rater of the
% disagreed FOG
FOG_summed_v2=FOG_vector{1}+2*FOG_vector{2}; % 0=agreed no FOG; 3=agreed FOG; 1=FOG only annotated by rater 1; 2=FOG only annotated by rater 2

% make an agreement table for this file
varnames={'subject', 'filename', 'nrFOG_rater1', 'durFOG_rater1', 'nrFOG_rater2', 'durFOG_rater2', 'durFOG_agreed', 'durFOG_disagreed_rater1', 'durFOG_disagreed_rater2', 'total_duration', 'kappa'};
vartypes=[repmat({'string'}, [1,2]), repmat({'double'}, [1,9])];

agreement_info.subject={ID};
[path, name, ext]=fileparts(filename_combined);
agreement_info.filename={name};
agreement_info.nrFOG_rater1=height(FOG_annotations{1});
agreement_info.durFOG_rater1=sum([FOG_annotations{1}.EndTime_Ss_msec-FOG_annotations{1}.BeginTime_Ss_msec]);
agreement_info.nrFOG_rater2=height(FOG_annotations{2});
agreement_info.durFOG_rater2=sum([FOG_annotations{2}.EndTime_Ss_msec-FOG_annotations{2}.BeginTime_Ss_msec]);
agreement_info.durFOG_agreed=sum(FOG_summed_v2==3)/sf;
agreement_info.durFOG_disagreed_rater1=sum(FOG_summed_v2==1)/sf;
agreement_info.durFOG_disagreed_rater2=sum(FOG_summed_v2==2)/sf;
agreement_info.total_duration=duration_gait_tasks{1};

% calculate kappa correlation coefficient of this file
kappa=kappacoefficient(agreement_info);
agreement_info.kappa=kappa;

% display
fprintf('Agreement info of this file: \n')
display(table(varnames(3:end)', agreement_info{:,3:end}', 'VariableNames', {'annotation_info', 'value'}))

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
    
% OVERLAPPINGEVT
function [idx] = overlappingevt(annotations, beginsample, endsample)
% find the indices of annotation events that fall within the event with the
% given [beginsample endsample].
  idx=find(([annotations.BeginTime_sample]<=beginsample & [annotations.EndTime_sample]>beginsample) |... % annotation includes the beginsample
    ([annotations.BeginTime_sample]<endsample & [annotations.EndTime_sample]>=endsample) | ... % annotation includes the endsample
    ([annotations.BeginTime_sample]>=beginsample & [annotations.EndTime_sample]<endsample)); % annotation falls within the event

% KAPPACOEFFICIENT
function kappa=kappacoefficient(agreement_t)
total_duration=sum(agreement_t.total_duration);

a=sum(agreement_t.durFOG_agreed);
b=sum(agreement_t.durFOG_disagreed_rater1);
c=sum(agreement_t.durFOG_disagreed_rater2);
d=total_duration-a-b-c;

Po=(a+d)/total_duration;
Pyes=((a+c)/total_duration)*((a+b)/total_duration);
Pno=((b+d)/total_duration)*((c+d)/total_duration);

kappa=(Po-(Pyes+Pno))/(1-(Pyes+Pno));


% SPEARMANCORRELATION
function [corr_nrFOG, corr_durFOG]=spearmancorrelation(agreement_t)
% Number of FOG
corr_nrFOG = corr(agreement_t.nrFOG_rater1, agreement_t.nrFOG_rater2, 'Type', 'Spearman');
% plot
figure;
scatter(agreement_t.nrFOG_rater1, agreement_t.nrFOG_rater2, 'filled');
title(sprintf('Spearman correlation for number of FOG: %.02f', corr_nrFOG)); 
xlabel('rater1'), ylabel('rater2');
grid on; axis square; 
% FIXME: scaling

% Duration of FOG
corr_durFOG = corr(agreement_t.durFOG_rater1, agreement_t.durFOG_rater2, 'Type', 'Spearman');
% plot
figure; 
scatter(agreement_t.durFOG_rater1, agreement_t.durFOG_rater2, 'filled');
title(sprintf('Spearman correlation for FOG duration (sec): %.02f', corr_durFOG));
xlabel('rater1'), ylabel('rater2');
grid on; axis square; 
% FIXME: scaling


