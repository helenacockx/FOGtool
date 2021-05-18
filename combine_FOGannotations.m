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

% TO DISCUSS:
% - what if 'export multiple files as'
% - file naming 
% - what if FOG_Trigger and FOG_Type are not dependent tiers?/are not both
% in the same line
% - what do we use as the length of the video file (because this has
% implications for the calculation of kappa and Spearman's correlation)
% - decide about corrections... Have a setting for this? (e.g. strict vs
% broad)
% - option to visualize/always visualize? use time vector that corresponds
% to ELAN? but also depends on offset master media file (so uncheck)
% - is the extra tier 'check_type' and 'check_trigger' handy, or remove
% this whole column?
% - still needs to calculate the kappa correlation and Spearman
% correlation.

%% set-up:
clear all; close all;

folder_rater1='\\dcn-srv.science.ru.nl\dcn\biophysics\prompt\freezing_fnirs\data\processed\annotations\Helena';
folder_rater2='\\dcn-srv.science.ru.nl\dcn\biophysics\prompt\freezing_fnirs\data\processed\annotations\Yuli';
folder_combined='\\dcn-srv.science.ru.nl\dcn\biophysics\prompt\freezing_fnirs\data\processed\annotations\combined';

subjects={'PD11'};

sf=60; % choose sampling frequency
ts=(1/sf); % time steps 

tolerance_sec=2;    % Tolerance in sec
tolerance=tolerance_sec*sf;

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
  
  %%
  for f=1:length(files{1})
    %% collect the annotations of these files
    for i=1:2
      annotations=readtable(fullfile(files{i}(f).folder, files{i}(f).name), 'ReadVariableNames', 1, 'HeaderLines', 0);
      % only select FOG annotations 
      FOG_annotations{i}=annotations(~ismissing(annotations.FOG_Trigger(:))|~ismissing(annotations.FOG_Trigger(:)),:);
      NOTES{i}=annotations(~ismissing(annotations.NOTES(:)) & (ismissing(annotations.FOG_Trigger(:))&ismissing(annotations.FOG_Type(:))),:); % extra notes that were not annotated during a FOG
      % check if each FOG has been labeled with a FOG_Trigger and a FOG_Type
      if any(ismissing(FOG_annotations{i}.FOG_Trigger(:))|ismissing(FOG_annotations{i}.FOG_Type))
        warning('Not all FOG events were both annotated for FOG_Trigger and FOG_Type')
      end
    end
    
    %% convert annotations to boolean vectors based on the given sampling frequency
    endtime=max([FOG_annotations{1}.EndTime_Ss_msec; FOG_annotations{2}.EndTime_Ss_msec]); % FIXME: or use length of video file instead?
    t=[1:1/sf:endtime+5*sf]; % time vector (add 5 extra seconds to make sure that the boolvec_FOG goes back to zero after the last FOG)
    boolvec=zeros(1,length(t)); % boolean vector
    for i=1:2
      % add extra column with begin sample and end sample of FOG annotation
      FOG_annotations{i}.BeginTime_sample=round(FOG_annotations{i}.BeginTime_Ss_msec*sf+1);
      FOG_annotations{i}.EndTime_sample=round(FOG_annotations{i}.EndTime_Ss_msec*sf+1);
      % create a boolean vector with the FOG annotations of this rater
      boolvec_FOG=boolvec;
      for k=1:height(FOG_annotations{i})
        boolvec_FOG(FOG_annotations{i}.BeginTime_sample(k):FOG_annotations{i}.EndTime_sample(k))=1;
      end
      FOG_vector{i}=boolvec_FOG;
    end
    
    %% combine annotations
    FOG_summed=FOG_vector{1}+FOG_vector{2}; % 0=agreed no FOG; 2=agreed FOG; 1=disagreed
    
    % find begin and end samples of disagreed events
    beginsample=find(FOG_summed==1 & diff([0 FOG_summed])~=0);% find beginsamples when going from agreed to disagreed 
    endsample=find(FOG_summed==1 &  diff([FOG_summed 0])~=0); % find endsamples when going back from disagreed to agreed
    if length(beginsample)~= length(endsample)
      error
    end
    
    % correct disagreements based on tolerance window length (e.g. 2 sec)
    % A; Als bij een plateau van 1 zowel voor als na een 2 zit
    % plateau<=tolerantie, dan concensus = 2.
    % B; als bij een plateau van 1 zowel voor als na een 0 zit, en
    % plateau<=tolerantie, dan consensus = 0. ==> or 1?
    % C; annoteerder A annoteert het begin eerder dan annoteerder B
    % Als voor een plateau een 0 zit en daarna een 2, plateau<=.5*tolerantie,
    % dan consensus = 2.
    % C; annoteerder A annoteert het einde eerder dan annoteerder B
    % Als voor een plateau een 2 zit en daarna een 0, plateau<=.5*tolerantie,
    % dan consensus = 2. 
    FOG_corrected=FOG_summed;
    for k=1:length(beginsample)
      if FOG_summed(beginsample(k)-1)== 2 & FOG_summed(endsample(k)+1)==2 
        if beginsample(k)-endsample(k)<= tolerance
          FOG_corrected(beginsample(k):endsample(k))=2;
        else
          FOG_corrected(beginsample(k):endsample(k))=1;% or 0
        end
      elseif FOG_summed(beginsample(k)-1)== 0 & FOG_summed(endsample(k)+1)==0
        if endsample(k)-beginsample(k)<= tolerance
          FOG_corrected(beginsample(k):endsample(k))=1; % not agreed annoation (or 0)
        else
          FOG_corrected(beginsample(k):endsample(k))=1; % not agreed annotation
        end
      elseif FOG_corrected(beginsample(k)-1)== 0 & FOG_summed(endsample(k)+1)==2
        if endsample(k)-beginsample(k)<= tolerance %0.5*
          FOG_corrected(beginsample(k):endsample(k))=2;
        else
          FOG_corrected(beginsample(k):endsample(k))=1; % not agreed annotation
        end          
      elseif FOG_corrected(beginsample(k)-1)== 2 & FOG_summed(endsample(k)+1)==0
        if endsample(k)-beginsample(k)<= tolerance %0.5*
          FOG_corrected(beginsample(k):endsample(k))=2;
        else
          FOG_corrected(beginsample(k):endsample(k))=1; % not agreed annotation
        end
      end
    end
    
    % find the agreed and disagreed FOG's
    FOG_agreed=(FOG_corrected==2);
    FOG_disagreed=(FOG_corrected==1);
    
    %% visualize
    figure;
    ax(1)=subplot(4,1,1); plot(t, FOG_vector{1}, 'Color', '#EDB120'); ylim([-1 2]); ylabel('FOG rater 1'); title(files{1}(f).name)
    ax(2)=subplot(4,1,2); plot(t, FOG_vector{2},  'Color','#0072BD'); ylim([-1 2]);ylabel('FOG rater 2'); title(files{2}(f).name)
    ax(3)=subplot(4,1,3); plot(t, FOG_agreed,  'Color','#77AC30'); ylim([-1 2]);ylabel('FOG agreed')
    ax(4)=subplot(4,1,4); plot(t, FOG_disagreed, 'Color', '#A2142F'); ylim([-1 2]);ylabel('FOG disagreed')
    linkaxes(ax)
    xlabel('time (in sec)');
    
    %% convert FOG disagreed to a table and add the rater, trigger and type for this FOG
    [beginsample, endsample]=vec2event(FOG_disagreed);% find beginsample and endsample of each event
    n=length(beginsample);
    FOG_disagreed_t=table(beginsample', endsample',nan(n,1), cell(n,1), cell(n,1),cell(n,1), cell(n,1), cell(n,1), cell(n,1),cell(n,1), cell(n,1), 'VariableNames', {'BeginTime_sample', 'EndTime_sample','rater', 'FOG_agreed_Trigger', 'FOG_agreed_Type', 'FOG_disagreed_Trigger', 'FOG_disagreed_Type', 'check_trigger', 'check_type','NOTES_rater1', 'NOTES_rater2'});% convert into table
    for k=1:height(FOG_disagreed_t)
      % find overlapping events of the two raters (FIXME: make helper
      % function?)
      idx_rater1=find(([FOG_annotations{1}.BeginTime_sample]<=beginsample(k)&[FOG_annotations{1}.EndTime_sample]>beginsample(k))|([FOG_annotations{1}.BeginTime_sample]<endsample(k)&[FOG_annotations{1}.EndTime_sample]>=endsample(k)));
      idx_rater2=find(([FOG_annotations{2}.BeginTime_sample]<=beginsample(k)&[FOG_annotations{2}.EndTime_sample]>beginsample(k))|([FOG_annotations{2}.BeginTime_sample]<endsample(k)&[FOG_annotations{2}.EndTime_sample]>=endsample(k)));
      if length(idx_rater1)==1 & isempty(idx_rater2)
        FOG_disagreed_t.rater(k)=1;
        FOG_disagreed_t.FOG_disagreed_Trigger(k)=FOG_annotations{1}.FOG_Trigger(idx_rater1);
        FOG_disagreed_t.FOG_disagreed_Type(k)=FOG_annotations{1}.FOG_Type(idx_rater1);
        FOG_disagreed_t.NOTES_rater1(k)=FOG_annotations{1}.NOTES(idx_rater1);
      elseif isempty(idx_rater1) & length(idx_rater2)==1
        FOG_disagreed_t.rater(k)=2;
        FOG_disagreed_t.FOG_disagreed_Trigger(k)=FOG_annotations{2}.FOG_Trigger(idx_rater2);
        FOG_disagreed_t.FOG_disagreed_Type(k)=FOG_annotations{2}.FOG_Type(idx_rater2);   
        FOG_disagreed_t.NOTES_rater2(k)=FOG_annotations{2}.NOTES(idx_rater2);
      else
        error % the disagreed annotations should only be rated by one person and no the event should only overlap with one other event
      end
    end

    %% convert FOG agreed to a table and check whether trigger and type for this FOG of both raters was the same (if not combine both values and check_trigger/check_type=true)
    [beginsample, endsample]=vec2event(FOG_agreed);% find beginsample and endsample of each event
    n=length(beginsample);
    FOG_agreed_t=table(beginsample', endsample', nan(n,1), cell(n,1), cell(n,1), cell(n,1), cell(n,1), cell(n,1), cell(n,1),cell(n,1), cell(n,1), 'VariableNames', {'BeginTime_sample', 'EndTime_sample', 'rater', 'FOG_agreed_Trigger', 'FOG_agreed_Type', 'FOG_disagreed_Trigger', 'FOG_disagreed_Type',  'check_trigger', 'check_type', 'NOTES_rater1', 'NOTES_rater2'});% convert into table
    for k=1:height(FOG_agreed_t)
      % find overlapping events of the other raters
      idx_rater1=find(([FOG_annotations{1}.BeginTime_sample]<=beginsample(k)&[FOG_annotations{1}.EndTime_sample]>beginsample(k))|([FOG_annotations{1}.BeginTime_sample]<endsample(k)&[FOG_annotations{1}.EndTime_sample]>=endsample(k)));
      idx_rater2=find(([FOG_annotations{2}.BeginTime_sample]<=beginsample(k)&[FOG_annotations{2}.EndTime_sample]>beginsample(k))|([FOG_annotations{2}.BeginTime_sample]<endsample(k)&[FOG_annotations{2}.EndTime_sample]>=endsample(k)));
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
    FOG_all_t=[FOG_agreed_t; FOG_disagreed_t];
    % add timing in seconds (instead of samples)
    FOG_all_t.BeginTime_Ss_msec=(FOG_all_t.BeginTime_sample(:)-1)/sf;
    FOG_all_t.EndTime_Ss_msec=(FOG_all_t.EndTime_sample(:)-1)/sf;
    % add extra notes
    for k=1:height(NOTES{1})
      FOG_all_t.BeginTime_Ss_msec(end+1)=NOTES{1}.BeginTime_Ss_msec(k);
      FOG_all_t.EndTime_Ss_msec(end)=NOTES{1}.EndTime_Ss_msec(k);
      FOG_all_t.NOTES_rater1(end)=NOTES{1}.NOTES(k);
    end
    for k=1:height(NOTES{2})
      FOG_all_t.BeginTime_Ss_msec(end+1)=NOTES{2}.BeginTime_Ss_msec(k);
      FOG_all_t.EndTime_Ss_msec(end)=NOTES{2}.EndTime_Ss_msec(k);
      FOG_all_t.NOTES_rater2(end)=NOTES{2}.NOTES(k);
    end
    %% export
    writetable(FOG_all_t, fullfile(folder_combined, sprintf('sub-%s_file-%02d.tsv', subjects{s}, f)),  'FileType', 'text', 'Delimiter', '\t')
    
  end
end


%% HELPER FUNCTIONS
function     [beginsample, endsample]=vec2event(boolvec)
    tmp = diff([0 boolvec 0]);
    beginsample = find(tmp==+1);
    endsample = find(tmp==-1) - 1;
end