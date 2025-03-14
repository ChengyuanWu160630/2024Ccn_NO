function[output_uncertainty,output_alpha_pwe,output_f,output_alpha,output_pwe,congruency,EEG_deleted,high_indices_f,low_indices_f,high_indices_alpha,low_indices_alpha,high_indices_pwe,low_indices_pwe]=EEG_aligh_behavior(f,alpha,congruency,uncertainty,EEG,EEGdel,behavior_del,residue) 
for p=1:length(congruency)
    if congruency(p)==0
        congruency(p)=-1;
    else
        congruency(p)=1;
    end
end
tem_congruency=congruency;
exclude=80:80:640;
behavior_del(ismember(behavior_del,exclude))=[];

pwe=1:length(f);
f_NO=1:length(f);
tem_f=f;
tem_alpha=alpha;
tem_uncertainty=uncertainty;

tem_f(behavior_del+1)=[];
tem_alpha(behavior_del+1)=[];
tem_uncertainty(behavior_del+1)=[];

zscore_f=zscore(tem_f);
zscore_alpha=zscore(tem_alpha);
zscore_uncertainty=zscore(tem_uncertainty);

tem_congruency(behavior_del+1)=[];
zscore_pwe=tem_congruency.*zscore_f;
need_to_changeZ=~ismember(f_NO,behavior_del+1);
f(need_to_changeZ)=zscore_f; 
alpha(need_to_changeZ)=zscore_alpha;
pwe(need_to_changeZ)=zscore_pwe;
uncertainty(need_to_changeZ)=zscore_uncertainty;
behavior_del=union(behavior_del,behavior_del+1);
behavior_del(ismember(behavior_del,exclude))=[];

%% calculate residual target model variable
if residue==1
    [residuals]=calculate_residue(pwe,alpha,f,congruency,uncertainty,need_to_changeZ);
    f=residuals.f;
    alpha=residuals.alpha;
    pwe=residuals.pwe;
    alpha_pwe=residuals.alpha_pwe;
    uncertainty=residuals.uncertainty;
    median_value_f = median(f(residuals.input_num+1));
    median_value_alpha = median(alpha(residuals.input_num+1));
    median_value_pwe = median(pwe(residuals.input_num));
else
    median_value_f = median(zscore_f);
    median_value_alpha = median(zscore_alpha);
    median_value_pwe = median(zscore_pwe);
    alpha_pwe=1:length(f_NO);
    uncertainty=1:length(f_NO);
end

%% remove bad epoch and label trial
trialNO=1:640;
EEG_invalidNO=80:80:640;
delete80=ismember(trialNO,EEG_invalidNO);
for i=1:length(EEGdel)
    trialNO(EEGdel{i})=[];
    delete80(EEGdel{i})=[];
end
invalid_indices = ismember(trialNO, EEG_invalidNO);
trialNO(invalid_indices) = [];
EEG_delete80=find(delete80==1);
EEG_remove80 = pop_rejepoch( EEG, EEG_delete80 ,0);

behavior2eeg_deletedtrial=ismember(trialNO, (behavior_del));
trialNO(behavior2eeg_deletedtrial) = [];
corre_trialNO=trialNO+1;
behavior2eeg_deletedtrial=find(behavior2eeg_deletedtrial==1);
EEG_deleted = pop_rejepoch( EEG_remove80, behavior2eeg_deletedtrial ,0);
congruency=congruency(trialNO);

output_uncertainty=uncertainty(corre_trialNO);
output_f=f(corre_trialNO);
output_alpha=alpha(corre_trialNO);
output_alpha_pwe=alpha_pwe(corre_trialNO);
output_pwe=pwe(trialNO);
high_indices_f = find(output_f > median_value_f);
low_indices_f = find(output_f < median_value_f);
high_indices_alpha = find(output_alpha > median_value_alpha);
low_indices_alpha= find(output_alpha < median_value_alpha);
high_indices_pwe = find(output_pwe > median_value_pwe);
low_indices_pwe= find(output_pwe < median_value_pwe);
end