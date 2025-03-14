function[struct_theta,struct_ALPHA,struct_beta1,struct_beta2,A,B,C,D,E]=TF_LME_sti(subject,behavior_path,EEG_path)
load(fullfile(behavior_path,[subject,'.mat']));
del_list=[];
for i=1:length(RTs)
    del=find(RTs{i}==-1)+80*(i-1);
    del_list=[del_list,del];
end
congruency=[c{1:8}];
for i=1:8
    load(fullfile(behavior_path,subject,['sim',num2str(i),'.mat']))
    alpha((80*i-79):(80*i))=vEstimate;
    f((80*i-79):(80*i))=pEstimate;
    uncertainty((80*i-79):(80*i))=stdEstimate;
end
fn=[subject,'.set'];
data_path=EEG_path;
EEG = pop_loadset('filename',fn,'filepath',data_path);
firstdel=EEG.manual_refuse;
seconddel=EEG.auto_refuse;
EEGdel={firstdel,seconddel};
%% remove epoch
EEG_R=EEG;
f_R=f;
alpha_R=alpha;
congruency_R=congruency;
[~,~,f,alpha,pwe,congruency,EEG,~,~,~,~,~,~]=EEG_aligh_behavior(f,alpha,congruency,uncertainty,EEG,EEGdel,del_list,0);%gai
[R_uncertainty,R_alpha_pwe,R_f,R_alpha,R_pwe,~,~,~,~,~,~,~,~]=EEG_aligh_behavior(f_R,alpha_R,congruency_R,uncertainty,EEG_R,EEGdel,del_list,1);
EEG = pop_epoch( EEG, {  '2   '  }, [-1           3.002], 'newname', 'EEProbe continuous data epochs', 'epochinfo', 'yes');
[~, idx] = sort([EEG.chanlocs.X], 'descend');
EEG.chanlocs = EEG.chanlocs(idx);
EEG.data = EEG.data(idx,:,:);
struct_theta.chanlocs = EEG.chanlocs;
struct_ALPHA.chanlocs = EEG.chanlocs;
struct_beta1.chanlocs = EEG.chanlocs;
struct_beta2.chanlocs = EEG.chanlocs;
%% time-frequency-analysis
ft_data = eeglab2fieldtrip(EEG, 'preprocessing', 'none');
cfg = [];
cfg.channel = 'all';
cfg.method = 'tfr';
cfg.output = 'pow';
cfg.foi = 3:0.5:30;
cfg.toi = [-0.64:0.002:-0.36,0:0.002:1];
cfg.width = 3;
cfg.keeptrials  = 'yes';
subj_data_tf = ft_freqanalysis(cfg,ft_data);

cfg_baseline = [];
cfg_baseline.baseline = [-0.6 -0.4];
cfg_baseline.baselinetype = 'absolute'; 
subj_data_tf_blc = ft_freqbaseline(cfg_baseline, subj_data_tf);
clear subj_data_tf

cfg_sel = [];
cfg_sel.latency = [0 1];   

cfg_alpha = cfg_sel;
cfg_alpha.frequency = [8 12];
alpha_data = ft_selectdata(cfg_alpha, subj_data_tf_blc);

cfg_theta = cfg_sel;
cfg_theta.frequency = [4 7];
theta_data = ft_selectdata(cfg_theta, subj_data_tf_blc);

cfg_beta1 = cfg_sel;
cfg_beta1.frequency = [14 20];
beta1_data = ft_selectdata(cfg_beta1, subj_data_tf_blc);
 
cfg_beta2 = cfg_sel;
cfg_beta2.frequency = [21 30];
beta2_data = ft_selectdata(cfg_beta2, subj_data_tf_blc);
clear subj_data_tf_blc

cfg_avg = [];
cfg_avg.avgoverfreq = 'yes'; 

theta_data = ft_selectdata(cfg_avg, theta_data);
alpha_data = ft_selectdata(cfg_avg, alpha_data);
beta1_data = ft_selectdata(cfg_avg, beta1_data);
beta2_data = ft_selectdata(cfg_avg, beta2_data);

theta = permute(squeeze(theta_data.powspctrm), [2, 3, 1]); 
ALPHA = permute(squeeze(alpha_data.powspctrm), [2, 3, 1]); 
beta_1 = permute(squeeze(beta1_data.powspctrm), [2, 3, 1]);
beta_2 = permute(squeeze(beta2_data.powspctrm), [2, 3, 1]); 
clear theta_data alpha_data beta1_data beta2_data

ALPHA=zscore(ALPHA,1,3);
beta_2=zscore(beta_2,1,3);
theta=zscore(theta,1,3);
beta_1=zscore(beta_1,1,3);

struct_theta.data = theta;
struct_ALPHA.data = ALPHA;
struct_beta1.data = beta_1;
struct_beta2.data = beta_2;
%% model variable matrix
A=[R_alpha;f;congruency;pwe;alpha.*f;alpha.*congruency;alpha.*pwe];
B=[alpha;R_f;congruency;pwe;alpha.*f;alpha.*congruency;alpha.*pwe];
C=[alpha;f;congruency;R_pwe;alpha.*f;alpha.*congruency;alpha.*pwe];
D=[alpha;f;congruency;pwe;alpha.*f;alpha.*congruency;R_alpha_pwe];
E=[alpha;f;congruency;pwe;alpha.*f;alpha.*congruency;alpha.*pwe;R_uncertainty];