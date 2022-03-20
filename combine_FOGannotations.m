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
% - create new ELAN file; make sure to add the same video offsets
% - import CSV/tab-delimited text file: tick 'FOG_agreed_Trigger' (Annotation),
% 'FOG_agreed_Type' (Annotation), 'FOG_disagreed_Trigger' (Annotation),
% 'FOG_disagreed_Type' (Annotation), 'check_trigger' (Annotation),
% 'check_type' (Annotation), 'NOTES_rater1' (Annotation), 'NOTES_rater2'
% (Annotation), 'begintime_Ss_msec' (Begin Time), 'endtime_Ss_msec'
% (End time); specify first row of data =2; specify delimiter = tab; Skip
% empty cells, don't create empty annotations = true.
% 4. the third rater decides on all the annotations in FOG_disagreed_... to
% keep it or delete it. When an annotation is flagged with 'check_type' or
% 'check_trigger', the raters did not agree on the value of the annoation,
% and the value of the type or trigger should be checked respectively
% 5. the new file can be saved and exported again to a fileformat of
% choice. By combining the agreed and disagreed FOG's (which have a
% consensus now), you have all the FOG's.

% TO DO:
% - what if 'export multiple files as'
% - option to visualize/always visualize? use time vector that corresponds
% to ELAN? but also depends on offset master media file (so uncheck)
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
%     % check if each FOG has been labeled with a FOG_Trigger and a FOG_Type
%     idx=find(ismissing(FOG_annotations{i}.fog_trigger(:))|ismissing(FOG_annotations{i}.fog_type));
%     if ~isempty(idx)
%       if length(idx)==height(FOG_annotations{i}) % only fog trigger/type was annotated by this rater
%         if 
%       else
%         warning('Not all FOG events were both annotated for FOG_Trigger and FOG_Type for rater %.0d', i)
%         display(FOG_annotations{i}(idx,:))
%       end
%     end
    
    % calculate total duration based on the gait tasks
    gait_tasks{i}=annotations{i}(~ismissing(annotations{i}.gait_task(:)),:);
    if isempty(gait_tasks{i})
        warning('No gait_tasks were found for rater %.0d', i)
        duration_gait_tasks{i}=nan;
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
    warning('Only annotations of rater %.0d contained gait_tasks. Using those to calculate total duration', rater)
    total_duration = duration_gait_tasks{rater};
    endtime = max(gait_tasks{rater}.endtime_msec);
    gait_tasks=gait_tasks{rater};
elseif isempty(gait_tasks{1}) & isempty(gait_tasks{2})
    warning('No gait tasks were found for any of the raters. Assuming that the gait_task ended after the last FOG event. The agreement parameters will not be calculated and the visualization might be wrong.')
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

% Use column with msec for samples
for i=1:height(gait_tasks)
    boolvec_task(gait_tasks.begintime_msec(i)+1:gait_tasks.endtime_msec(i)+1)=0;
end

for i=1:2
    % add extra column with begin sample and end sample of FOG annotation, nu
    % overbodig?
    % FOG_annotations{i}.begintime_msec=round(FOG_annotations{i}.begintime_Ss_msec*sf+1);
    % FOG_annotations{i}.endtime_msec=round(FOG_annotations{i}.endtime_Ss_msec*sf+1);
    
    % create a boolean vector with the FOG annotations of this rater
    boolvec_FOG=boolvec_task;
    for k=1:height(FOG_annotations{i})
        boolvec_FOG((FOG_annotations{i}.begintime_msec(k)+1):(FOG_annotations{i}.endtime_msec(k)+1))=1;
    end
    % check if all FOGs are falling insside the gait_task
    if sum(boolvec_FOG ==1 & isnan(boolvec_task))>0
      warning('%d FOG events are falling outside the gait_task and are removed from the list', numel(find(diff([boolvec_FOG ==1 & isnan(boolvec_task)])==1)))
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
varnames = {'begintime_msec', 'endtime_msec','rater', 'FOG_agreed_Trigger', ...
    'FOG_agreed_Type', 'FOG_disagreed_Trigger', 'FOG_disagreed_Type', 'check_trigger', 'check_type','NOTES_rater1', 'NOTES_rater2'};
vartypes(1,[1:3])={'double'};
vartypes(1,[4:11])={'string'};
FOG_disagreed_t=table('Size', [n, length(vartypes)], 'VariableNames', varnames, 'VariableTypes', vartypes);
FOG_disagreed_t.begintime_msec=beginsample'-1;
FOG_disagreed_t.endtime_msec=endsample'-1;

% find rater, trigger and type of this FOG
for k=1:height(FOG_disagreed_t)
    % find annotations that fall within this event
    idx_rater1 = overlappingevt(FOG_annotations{1}, beginsample(k), endsample(k));
    idx_rater2 = overlappingevt(FOG_annotations{2}, beginsample(k), endsample(k));
    
    if length(idx_rater1)==1 & isempty(idx_rater2) % this is an annotation from rater 1       
        FOG_disagreed_t.rater(k)=1;
        FOG_disagreed_t.FOG_disagreed_Trigger(k)=FOG_annotations{1}.fog_trigger(idx_rater1);
        FOG_disagreed_t.FOG_disagreed_Type(k)=FOG_annotations{1}.fog_type(idx_rater1);
        try
            FOG_disagreed_t.NOTES_rater1(k)=FOG_annotations{1}.notes(idx_rater1);
        end
    elseif isempty(idx_rater1) & length(idx_rater2)==1 % this is an annotation from rater 2
        FOG_disagreed_t.rater(k)=2;
        FOG_disagreed_t.FOG_disagreed_Trigger(k)=FOG_annotations{2}.fog_trigger(idx_rater2);
        FOG_disagreed_t.FOG_disagreed_Type(k)=FOG_annotations{2}.fog_type(idx_rater2);
        try
            FOG_disagreed_t.NOTES_rater2(k)=FOG_annotations{2}.notes(idx_rater2);
        end
    else
        error('Multiple annotations were found for this disagreed FOG episode') 
    end
end

%% convert FOG agreed to a table and check whether trigger and type for this FOG of both raters was the same (if not combine both values and check_trigger/check_type=true)
[beginsample, endsample]=vec2event(FOG_agreed);% find beginsample and endsample of each event
n=length(beginsample);
FOG_agreed_t=table('Size', [n, length(vartypes)], 'VariableNames', varnames, 'VariableTypes', vartypes);
FOG_agreed_t.begintime_msec=beginsample'-1;
FOG_agreed_t.endtime_msec=endsample'-1;

for k=1:height(FOG_agreed_t)
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
        FOG_agreed_t.FOG_agreed_Trigger(k)=triggers(~ismissing(triggers));
    else % a different value was given for trigger
        FOG_agreed_t.check_trigger(k)={'check_trigger'};
        % combine both values to one string
        trig_combi=triggers{1};
        for j=2:length(triggers)
            trig_combi=[trig_combi ' / ' triggers{j}];
        end
        FOG_agreed_t.FOG_agreed_Trigger(k)={trig_combi};
    end
    % check type
    types=[FOG_annotations{1}.fog_type(idx_rater1); FOG_annotations{2}.fog_type(idx_rater2)];
    types=unique(types);% unique values for types
    if length(types)==1 | flag_notype % the same value was given for FOG_Type/only one annotator characterized the type
        FOG_agreed_t.FOG_agreed_Type(k)=types(~ismissing(types));
    else  % a different value was given for type
        FOG_agreed_t.check_type(k)={'check_type'};
        % combine both values to one string
        type_combi=types{1};
        for h=2:length(types)
            type_combi=[type_combi ' / ' types{h}];
        end
        FOG_agreed_t.FOG_agreed_Type(k)={type_combi};
    end
    try
        FOG_agreed_t.NOTES_rater1(k)=FOG_annotations{1}.notes(idx_rater1);
    end
    try
        FOG_agreed_t.NOTES_rater2(k)=FOG_annotations{2}.notes(idx_rater2);
    end
end
%% Picture the results
color_discuss = sscanf('9b9b9b','%2x%2x%2x',[1 3])/255;
[~, name, ~]=fileparts(filename_combined);
setPos=[-7.5, 0.4, 0];    % position of y labels
% Turn warning about excluding colorbars, legends and non-axes off
warning('off','MATLAB:linkaxes:RequireDataAxes'); 
fntsz=20;

figure(4)
ax(1)=subplot(10,1,1);
for N=1:height(gait_tasks)
  clear t2 FilledArea timedata gait_stukje
        timedata=t(gait_tasks.begintime_msec(N):gait_tasks.endtime_msec(N));
        gait_stukje=ones(size(timedata));
        t2=[timedata, fliplr(timedata)];
        FilledArea=[gait_stukje, zeros(size(timedata))];
        fill(t2, FilledArea, 'k','LineStyle','none')
        hold on  
end
hold off
xlim([t(1) t(end)])
ylabel('Gait task','rotation', 0,'HorizontalAlignment','left',...
    'pos', [setPos(1), 0, setPos(2)]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'FontSize',fntsz);

ax(2)=subplot(11,1,[2,3,4]);
if length(FOG_annotations{1}.begintime_msec) == 0
    % Don't plot, no FOG annotated
else
    for M=1:length(FOG_annotations{1}.begintime_msec)       % plot ieder stukje afzonderlijk
        clear t2 FilledArea timedata FOG_stukje
        timedata=t(FOG_annotations{1}.begintime_msec(M):FOG_annotations{1}.endtime_msec(M));
        FOG_stukje=ones(size(timedata));
        t2=[timedata, fliplr(timedata)];
        FilledArea=[FOG_stukje, zeros(size(timedata))];
        fill(t2, FilledArea, 'k','LineStyle','none')
        hold on
    end
    hold off
end
ylabel('Rater 1','rotation', 0,'HorizontalAlignment','left',...
    'pos', setPos);
set(gca,'FontSize',fntsz);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
xlim([t(1) t(end)])
ylim ([0 1])

ax(3)=subplot(11,1,[5,6,7]);
if length(FOG_annotations{1}.begintime_msec) == 0
    % Don't plot, no FOG annotated
else
    for K=1:length(FOG_annotations{2}.begintime_msec)
        clear t2 FilledArea timedata FOG_stukje
        timedata=t(FOG_annotations{2}.begintime_msec(K):FOG_annotations{2}.endtime_msec(K));
        FOG_stukje=ones(size(timedata));
        t2=[timedata, fliplr(timedata)];
        FilledArea=[FOG_stukje, zeros(size(timedata))];
        fill(t2, FilledArea, 'k','LineStyle','none')
        hold on
    end
    hold off
end
ylabel('Rater 2','rotation', 0,'HorizontalAlignment','left',...
    'pos', setPos);
set(gca,'FontSize',fntsz);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
xlim([t(1) t(end)])
ylim ([0 1])

ax(4)=subplot(11,1,[8:10]);
%nexttile([3 1])
if length(FOG_agreed_t.begintime_msec) == 0
    % Don't plot, no FOG annotated
else
    for P=1:length(FOG_agreed_t.begintime_msec)
        clear t2 FilledArea timedata FOG_stukje
        timedata=t(FOG_agreed_t.begintime_msec(P):FOG_agreed_t.endtime_msec(P));
        FOG_stukje=ones(size(timedata));
        t2=[timedata, fliplr(timedata)];
        FilledArea=[FOG_stukje, zeros(size(timedata))];
        fill(t2, FilledArea, 'k','LineStyle','none')
        hold on
    end
    if length(FOG_disagreed_t.begintime_msec) == 0
        % Don't plot, no disagreed FOG episodes
    else
        for L=1:length(FOG_disagreed_t.begintime_msec)
            clear t2 FilledArea timedata FOG_stukje
            timedata=t(FOG_disagreed_t.begintime_msec(L):FOG_disagreed_t.endtime_msec(L));
            FOG_stukje=ones(size(timedata));
            t2=[timedata, fliplr(timedata)];
            FilledArea=[FOG_stukje, zeros(size(timedata))];
            
            %ax(3)=subplot(3,1,3);
            fill(t2, FilledArea, color_discuss,'LineStyle','none')
            hold on
        end
        
        hold off
    end
end
ylabel({'Ratings', 'combined'},'rotation', 0,'HorizontalAlignment','left',...
    'pos', [setPos(1), 0.25, setPos(3)]);
set(gca,'FontSize',fntsz);
set(gca,'YTickLabel',[]);
xlim([t(1) t(end)])
ylim ([0 1])
xlabel('Time (in seconds)');
linkaxes(ax);
%% combine the agreed and disagreed tables into one table
% add timing in seconds (instead of samples)?? Required?
% Check voor de -1
FOG_agreed_t.begintime_msec=(FOG_agreed_t.begintime_msec(:));
FOG_agreed_t.endtime_msec=(FOG_agreed_t.endtime_msec(:));
FOG_disagreed_t.begintime_msec=(FOG_disagreed_t.begintime_msec(:));
FOG_disagreed_t.endtime_msec=(FOG_disagreed_t.endtime_msec(:));
gait_tasks.begintime_msec=gait_tasks.begintime_msec;
gait_tasks.endtime_msec=gait_tasks.endtime_msec;

FOG_all_t=[FOG_agreed_t; FOG_disagreed_t];

%% Rearrange for easy import in ELAN
clear varnames vartypes
varnames = {'Tier','begintime_msec','endtime_msec', 'Annotation'};
vartypes(1,[2,3])={'double'};
vartypes(1,[1,4])={'string'};
FOG_Compared=table('Size', [1, length(vartypes)], 'VariableNames', varnames, 'VariableTypes', vartypes);
Extra_rows=table('Size', [1, length(vartypes)], 'VariableNames', varnames, 'VariableTypes', vartypes);

% Supress warning
warning('off','MATLAB:table:RowsAddedExistingVars');

% Add all freezing to the export
for KL=1:height(FOG_all_t)
    clear r
    [r,~]=size(FOG_Compared);
    if r == 1
        r=0;     % for the start of the table
    end
    % identify cells without data per row
    TF = ~ismissing(FOG_all_t(KL,:),{'' '.' [] NaN missing ""});
    indxTF=find(TF);
    indxTF=indxTF(3:end);
    datawaardes=FOG_all_t(KL,indxTF);
    FOG_Compared.Tier(r+1:r+width(datawaardes))=datawaardes.Properties.VariableNames';
    FOG_Compared.Annotation(r+1:r+width(datawaardes))=table2array(datawaardes(1,1:width(datawaardes)));
    FOG_Compared.begintime_msec(r+1:r+width(datawaardes),1)=ones(width(datawaardes),1).*table2array(FOG_all_t(KL,1));
    FOG_Compared.endtime_msec(r+1:r+width(datawaardes),1)=ones(width(datawaardes),1).*table2array(FOG_all_t(KL,2));
end

%% Add gait tasks
if ~flag_nogaittask
for LM=1:height(gait_tasks)
    clear r TF
    [r,~]=size(FOG_Compared);
    if r == 1
        r=0;     % for the start of the table
    end
    % identify cells without data per row
    TF = ~ismissing(gait_tasks(LM,:),{'' '.' [] NaN missing ""});
    indxTF=find(TF);
    indxTF=indxTF(3:end);
    datawaardes2=gait_tasks(LM,indxTF);
    FOG_Compared.Tier(r+1:r+width(datawaardes2))=datawaardes2.Properties.VariableNames';
    FOG_Compared.Annotation(r+1:r+width(datawaardes2))=table2array(datawaardes2(1,1:width(datawaardes2)));
    FOG_Compared.begintime_msec(r+1:r+width(datawaardes2),1)=ones(width(datawaardes2),1).*table2array(gait_tasks(LM,1));
    FOG_Compared.endtime_msec(r+1:r+width(datawaardes2),1)=ones(width(datawaardes2),1).*table2array(gait_tasks(LM,2));
end
end
% Check the headers for tiers in files besides previously stated to add to export
extra_rater1 = table_extra(annotations{1,1}, Extra_rows);
extra_rater2 = table_extra(annotations{1,2}, Extra_rows);

% Check for duplicate annotations
indx_raters=find(ismember(extra_rater1,extra_rater2));
extra_rater1(indx_raters,:)=[];
extra_raters = [extra_rater1; extra_rater2];

FOG_Compared=[FOG_Compared; extra_raters];

%% Write notes to file if notes are not empty.
NOTES_rater1 = tableNOTES(annotations{1,1}, Extra_rows, 1);
NOTES_rater2 = tableNOTES(annotations{1,2}, Extra_rows, 2);

if height(NOTES_rater1) == 1 && NOTES_rater1.endtime_msec == 0
    % notes are empty
else 
FOG_Compared=[FOG_Compared; NOTES_rater1];
end

if height(NOTES_rater2) == 1 && NOTES_rater2.endtime_msec == 0
    % notes are empty
else 
FOG_Compared=[FOG_Compared; NOTES_rater2];
end

writetable(FOG_Compared, filename_combined,  'FileType', 'text', 'Delimiter', '\t')

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
agreement_info.duration_FOG_rater1=sum([FOG_annotations{1}.endtime_msec-FOG_annotations{1}.begintime_msec]);
agreement_info.number_FOG_rater2=height(FOG_annotations{2});
agreement_info.duration_FOG_rater2=sum([FOG_annotations{2}.endtime_msec-FOG_annotations{2}.begintime_msec]);
agreement_info.number_FOG_agreed=height(FOG_agreed_t); % only for info, not to calculate agreement (because uses adjusted FOG annotations)
agreement_info.duration_FOG_agreed=sum(FOG_summed_v2==3)/sf;
agreement_info.number_FOG_disagreed_rater1=height(FOG_disagreed_t(FOG_disagreed_t.rater==1,:)); % only for info, not to calculate agreement (because uses adjusted FOG annotations)
agreement_info.duration_FOG_disagreed_rater1=sum(FOG_summed_v2==1)/sf;
agreement_info.number_FOG_disagreed_rater2=height(FOG_disagreed_t(FOG_disagreed_t.rater==2,:)); % only for info, not to calculate agreement (because uses adjusted FOG annotations)
agreement_info.duration_FOG_disagreed_rater2=sum(FOG_summed_v2==2)/sf;
agreement_info.total_duration=total_duration;

[agreement_info.positive_agreement, agreement_info.negative_agreement,...
    agreement_info.prevalence_index] = agreementParameters(agreement_info);
  
if flag_nogaittask
  warning('No gait tasks were found, do not calculate agreement parameters')
  agreement_info.total_duration = nan;
  agreement_info.positive_agreement = nan;
  agreement_info.negative_agreement = nan;
  agreement_info.prevalence_index = nan;
end

% calculate kappa correlation coefficient of this file
% kappa=kappacoefficient(agreement_info);
% agreement_info.kappa=kappa;

% calculate %agreement on trigger and type
if flag_notrigger
  agreement_info.agreement_trigger = nan;
else
  agreement_info.agreement_trigger = sum(~strcmp(FOG_agreed_t.check_trigger, 'check_trigger'))/height(FOG_agreed_t);
end
if flag_notype
  agreement_info.agreement_type = nan;
else 
  agreement_info.agreement_type = sum(~strcmp(FOG_agreed_t.check_type, 'check_type'))/height(FOG_agreed_t);
end

% display
fprintf('Agreement info of file %s: \n', name)
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

% TABLE_EXTRA
function Extra_rows = table_extra(annotation, Extra_rows)
% Check if other tiers are added to the annotations
listVarAnn1=annotation.Properties.VariableNames;
listVarKnown={'begintime_msec','endtime_msec','gait_task','FOG_Trigger','FOG_Type','NOTES'};
match1  = find(ismember(listVarAnn1, listVarKnown) == 0);
for G=1:length(match1)
    indxdata1=find(~ismissing(annotation.(listVarAnn1{match1(G)})));
    
    for H = 1:length(indxdata1)
        [r,~]=size(Extra_rows);
        if r == 1
            r = 0;     % for the start of the table
        end
        Extra_rows.Tier(r+H) = listVarAnn1{match1(G)};
        Extra_rows.Annotation(r+H) = annotation.(listVarAnn1{match1(G)})(indxdata1(H));
        Extra_rows.begintime_msec(r+H) = annotation.begintime_msec(indxdata1(H));
        Extra_rows.endtime_msec(r+H) = annotation.endtime_msec(indxdata1(H));
    end
end

% TABLE_NOTES
function Extra_rows = tableNOTES(annotation, Extra_rows, raternmbr)
indxdata1=find(~ismissing(annotation.NOTES));
if isempty(indxdata1)
    fprintf('%s%d\n','No notes were created by rater ', raternmbr)
else
for H = 1:length(indxdata1)
    [r,~]=size(Extra_rows);
    if r == 1
        r = 0;     % for the start of the table
    end
    Extra_rows.Tier(r+H) = {sprintf('%s%d', 'NOTES_rater', raternmbr)};
    Extra_rows.Annotation(r+H) = annotation.(listVarAnn1{match1})(indxdata1(H));
    Extra_rows.begintime_msec(r+H) = annotation.begintime_msec(indxdata1(H));
    Extra_rows.endtime_msec(r+H) = annotation.endtime_msec(indxdata1(H));
end
end


% OVERLAPPINGEVT
function [idx] = overlappingevt(annotations, beginsample, endsample)
% find the indices of annotation events that fall within the event with the
% given [beginsample endsample].
idx=find(([annotations.begintime_msec]<=beginsample & [annotations.endtime_msec]>beginsample) |... % annotation includes the beginsample
    ([annotations.begintime_msec]<endsample & [annotations.endtime_msec]>=endsample) | ... % annotation includes the endsample
    ([annotations.begintime_msec]>=beginsample & [annotations.endtime_msec]<endsample)); % annotation falls within the event

% AGREEMENTPARAMETERS
function [pos_agree, neg_agree,prev_indx] = agreementParameters(agreement_t)
n=sum(agreement_t.total_duration);
a=sum(agreement_t.duration_FOG_agreed);
b=sum(agreement_t.duration_FOG_disagreed_rater1);
c=sum(agreement_t.duration_FOG_disagreed_rater2);
d=n-a-b-c;

pos_agree =2*a/(n+(a-d));
neg_agree = 2*d/(n-(a-d));

prev_indx =(a-d)/n;



% KAPPACOEFFICIENT
function kappa =kappacoefficient(agreement_t)
total_duration=sum(agreement_t.total_duration);

a=sum(agreement_t.durFOG_agreed);
b=sum(agreement_t.durFOG_disagreed_rater1);
c=sum(agreement_t.durFOG_disagreed_rater2);
d=total_duration-a-b-c;

Po=(a+d)/total_duration;
Pyes=((a+c)/total_duration)*((a+b)/total_duration);
Pno=((b+d)/total_duration)*((c+d)/total_duration);

kappa=(Po-(Pyes+Pno))/(1-(Pyes+Pno));


