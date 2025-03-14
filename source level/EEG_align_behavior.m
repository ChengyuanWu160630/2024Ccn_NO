function[UN_f,behavior_uncertainty,alpha_pwe,behavior_f,behavior_alpha,behavior_pwe,congruency,EEG_deleted]=EEG_align_behavior(f,alpha,congruency,uncertainty,EEG,EEGdel,behavior_del,ridue) 
for p=1:length(congruency)
    if congruency(p)==0
        congruency(p)=-1;
    else
        congruency(p)=1;
    end
end
congruency2=congruency;
exclude=80:80:640;
behavior_del(ismember(behavior_del,exclude))=[];
pwe=1:length(f);
f_NO=1:length(f);
a_PE=1:length(f);
Z_f=f;
Z_alpha=alpha;
Z_uncertainty=uncertainty;
copyalpha=alpha;
copyF_U=f./uncertainty;
Z_F_U=copyF_U;
Z_f(behavior_del+1)=[];
Z_alpha(behavior_del+1)=[];
Z_uncertainty(behavior_del+1)=[];
Z_F_U(behavior_del+1)=[];
congruency2(behavior_del+1)=[];
Sec_PE=(Z_f-mean(Z_f))./0.5;
Sec_PE = max(Sec_PE, -1); 
Sec_PE= min(Sec_PE, 1); 
Sec_PE=1-(congruency2.*Sec_PE);
aSec_PE=Sec_PE;

zscore_f=zscore(Z_f);
zscore_alpha=zscore(Z_alpha);
zscore_uncertainty=zscore(Z_uncertainty);
zscore_F_U=zscore(Z_F_U);

zscore_pwe=congruency2.*zscore_f;
need_to_changeZ=~ismember(f_NO,behavior_del+1);
f(need_to_changeZ)=zscore_f; 
alpha(need_to_changeZ)=zscore_alpha;
pwe(need_to_changeZ)=zscore_pwe;
uncertainty(need_to_changeZ)=zscore_uncertainty;
copyF_U(need_to_changeZ)=zscore_F_U;
a_PE(need_to_changeZ)=aSec_PE;
behavior_del=union(behavior_del,behavior_del+1);
behavior_del(ismember(behavior_del,exclude))=[];
%% calculate residue
if ridue==1
    [residuals]=Scalculate_residue2(pwe,alpha,f,congruency,uncertainty,copyF_U,a_PE,copyalpha,need_to_changeZ);
    f=residuals.f;
    alpha=residuals.alpha;
    pwe=residuals.pwe;
    alpha_pwe=residuals.alpha_pwe;
    uncertainty=residuals.uncertainty;
    UN_f=residuals.uncertainty_f;
    median_value_f = median(f(residuals.input_num+1));
    median_value_alpha = median(alpha(residuals.input_num+1));
    median_value_pwe = median(pwe(residuals.input_num));
else
    median_value_f = median(zscore_f);
    median_value_alpha = median(zscore_alpha);
    median_value_pwe = median(zscore_pwe);
    alpha_pwe=1:length(f_NO);
    UN_f=1:length(f_NO);
end
%% 预测模式下剔除脑电分段和行为数据
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
behavior_uncertainty=uncertainty(corre_trialNO);
behavior_f=f(corre_trialNO);
behavior_alpha=alpha(corre_trialNO);
alpha_pwe=alpha_pwe(corre_trialNO);
UN_f=UN_f(corre_trialNO);
behavior_pwe=pwe(trialNO);
end
