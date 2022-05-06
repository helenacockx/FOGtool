
function agreement = agreement_calculator(agreement_t)
n=sum(agreement_t.total_duration);
a=sum(agreement_t.duration_FOG_agreed);
b=sum(agreement_t.duration_FOG_disagreed_rater1);
c=sum(agreement_t.duration_FOG_disagreed_rater2);
d=n-a-b-c;

% agreement parameters
agreement.pos_agree =2*a/(n+(a-d));
agreement.neg_agree = 2*d/(n-(a-d));
agreement.prev_indx =(a-d)/n;

% kappa
% Po=(a+d)/n;
% Pc=(((a+c)*(a+b))/n + ((b+d)*(c+d))/n)/n;
% agreement.kappa=(Po-Pc)/(1-Pc);
