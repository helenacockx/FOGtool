folder_combined='\\dcn-srv.science.ru.nl\dcn\biophysics\prompt\freezing_fnirs\data\processed\annotations\combined';
agreement_t_file=fullfile(folder_combined, 'agreement_table.tsv');

agreement_t = readtable(agreement_t_file, 'FileType', 'text');

% Kappa correlation
total_duration=sum(agreement_t.total_duration);

a=sum(agreement_t.durFOG_agreed);
b=sum(agreement_t.durFOG_disagreed_rater1);
c=sum(agreement_t.durFOG_disagreed_rater2);
d=total_duration-a-b-c;

Po=(a+d)/total_duration;
Pyes=((a+c)/total_duration)*((a+b)/total_duration);
Pno=((b+d)/total_duration)*((c+d)/total_duration);

kappa=(Po-(Pyes+Pno))/(1-(Pyes+Pno));
pabak= 2*Po - 1; % prevalence-ajusted bias-adjusted kappa

% intraclass-correlation coefficient
