function MVGC_coupling_beta(subject)
datapath=['/home/ChenQi_01/Downloads/ningbo_data/ningbo_behavior/behavior_result/',subject];
load(fullfile(datapath,[subject,'.mat']));
del_list=[];
for i=1:length(RTs)
    del=find(RTs{i}==-1)+80*(i-1);
    del_list=[del_list,del];
end
congruency=[c{1:8}];
for i=1:8
    load(strcat(datapath,'/sim',num2str(i),'.mat'))
    alpha((80*i-79):(80*i))=vEstimate;
    f((80*i-79):(80*i))=pEstimate;
    uncertainty((80*i-79):(80*i))=stdEstimate;
end
fn=[subject,'.set'];
data_path=['/home/ChenQi_01/Downloads/ningbo_data/finished'];
EEG = pop_loadset('filename',fn,'filepath',data_path);
firstdel=EEG.manual_refuse;
seconddel=EEG.auto_refuse;
EEGdel={firstdel,seconddel};
%% remove epoch
[~,~,~,~,~,~,~,EEG]=EEG_align_behavior(f,alpha,congruency,uncertainty,EEG,EEGdel,del_list,0);
EEG = pop_epoch( EEG, {  '2   '  }, [-1           3.002], 'newname', 'EEProbe continuous data epochs', 'epochinfo', 'yes');
%% 
% define significant electrodes
group1 = {'E32','E26','E28','E29','E33','E34','E35','E36','E39','E40','E45'};
group2 = {'E86','E92','E97','E98','E93','E87'};
group3 = {'E1','E2','E8','E9','E122','E123','E124','E117','E116','E110','E111','E109','E108','E104','E103'};
group4 = {'E47','E42','E52','E53','E60','E66','E54','E61','E67','E62','E72','E78','E77','E84','E85','E90','E91'};
group5 ={'E112','E106','E105','E80','E79','E6','E37','E16'};
group6 ={'E86','E87','E92','E93','E104','E105','E106','E110','E111','E112','E117','E118','E124'};
all_electrodes = [group1, group2, group3, group4,group5];
unique_electrodes= unique([group1, group2, group3, group4,group5,group6]);
EEG_beta2= pop_select( EEG, 'channel',group6);
EEG_all= pop_select( EEG, 'channel',unique_electrodes);
EEG = pop_select( EEG, 'channel',all_electrodes);
groups = {group1, group2, group3, group4, group5};
pairs = {};

for i = 1:length(groups)
    j=6;
        if j > i
         
            g1 = groups{i};
            g2 = group6;          
           
            for k = 1:length(g1)
                for l = 1:length(g2)
                    pairs{end+1, 1} = g1{k};
                    pairs{end, 2} = g2{l};
                end
            end
        end
    
end
pairs = reshape(pairs, [], 2);
labelslist1 = {EEG.chanlocs.labels};
labelslist2 = {EEG_beta2.chanlocs.labels};
index_pairs = zeros(size(pairs));
for i = 1:size(pairs, 1)
    for j = 1:size(pairs, 2)
        if j==1
            index_pairs(i, j) = find(strcmp(labelslist1, pairs{i, j}));
        elseif j==2
            index_pairs(i, j) = find(strcmp(labelslist2, pairs{i, j}));
        end
    end
end

%% substitute data
ft_data1 = eeglab2fieldtrip(EEG, 'preprocessing', 'none');
ft_data2 = eeglab2fieldtrip(EEG_beta2, 'preprocessing', 'none');
ft_data_all = eeglab2fieldtrip(EEG_all, 'preprocessing', 'none');
load(fullfile('/home/ChenQi_01/Downloads/method_wavelet/direction_connection/granger/nonpara/coupling/FRE_nopading',subject,'FRE_3band(beta).mat'));
ALPHA_SPM=squeeze(alpha_data.powspctrm);
theta_SPM=squeeze(theta_data.powspctrm);
beta2_SPM=squeeze(beta2_data.powspctrm);
ALPHA_SPM(isnan(ALPHA_SPM))=0;
theta_SPM(isnan(theta_SPM))=0;
beta2_SPM(isnan(beta2_SPM))=0;
for i=1:length(ft_data1.trial)
    for j=1:length(ft_data1.label)
        index=find(strcmp(ft_data_all.label,ft_data1.label{j}));
        if sum(strcmp(group5,ft_data1.label{j}))>0
            ft_data1.trial{i}(j,:)=squeeze(theta_SPM(i,index,:));
        else
            ft_data1.trial{i}(j,:)=squeeze(ALPHA_SPM(i,index,:));
        end
    end
end
for i=1:length(ft_data2.trial)
    for j=1:length(ft_data2.label)
        index=find(strcmp(ft_data_all.label,ft_data2.label{j}));
        ft_data2.trial{i}(j,:)=squeeze(beta2_SPM(i,index,:));
    end
end
cfg=[];
cfg.latency=[-0.63 2.63];
ft_data1=ft_selectdata(cfg,ft_data1);
ft_data2=ft_selectdata(cfg,ft_data2);
cfg=[];
cfg.avgoverrpt ='yes';
ft_data1_avg=ft_selectdata(cfg,ft_data1);
ft_data2_avg=ft_selectdata(cfg,ft_data2);
EEG = pop_epoch( EEG, {  '2   '  }, [-0.63           2.632], 'newname', 'EEProbe continuous data epochs', 'epochinfo', 'yes');
EEG_beta2 = pop_epoch( EEG_beta2, {  '2   '  }, [-0.63           2.632], 'newname', 'EEProbe continuous data epochs', 'epochinfo', 'yes');
for i=1:length(ft_data1.trial)
    EEG.data(:,:,i)=ft_data1.trial{i}-ft_data1_avg.trial{1};
    EEG_beta2.data(:,:,i)=ft_data2.trial{i}-ft_data2_avg.trial{1};
%     EEG.data(:,:,i)=ft_data1.trial{i};
%     EEG_beta2.data(:,:,i)=ft_data2.trial{i};
end
EEG = pop_resample( EEG, 200);
EEG_beta2 = pop_resample( EEG_beta2, 200);
%% Parameters
regmode   = 'LWR';   
fs        = 200;    
fres      = 1000;    
morder    = 15;      
win= 100;
tstat     = 'F';   
%% Vertical regression GC calculation
glanger_time=round([-0.35:0.01:1]*1000);
timedomain_G=zeros([2,length(glanger_time),size(pairs,1)],'single');
for elec=1:size(pairs,1)
    idx=index_pairs(elec,:);
    X=zeros([2,size(squeeze(EEG.data(idx(1),:,:)))]);
    X(1,:,:)=squeeze(EEG.data(idx(1),:,:));
    X(2,:,:)=squeeze(EEG_beta2.data(idx(2),:,:));
    for e = 1:length(glanger_time)
        j = find(EEG.times==glanger_time(e));
        slidewin=X(:,j-win/2:j+win/2,:);

        for T=1:size(slidewin,3)
            slidewin(1,:,T)=zscore(detrend(squeeze(slidewin(1,:,T))));
            slidewin(2,:,T)=zscore(detrend(squeeze(slidewin(2,:,T))));
%           slidewin(1,:,T)=zscore(squeeze(slidewin(1,:,T)));
%           slidewin(2,:,T)=zscore(squeeze(slidewin(2,:,T)));
        end

    try
        [A,SIG] = tsdata_to_var(slidewin,morder,regmode);
    catch
        fprintf('Skipping time point %d due to an error.\n', e);
        continue;
    end
        info = var_info(A,SIG);
        [F,~] = var_to_pwcgc(A,SIG,slidewin,regmode,tstat);
    timedomain_G(1,e,elec)=F(2,1);
    timedomain_G(2,e,elec)=F(1,2);
    end
end
results_folder='/home/ChenQi_01/Downloads/method_wavelet/direction_connection/granger/mvgc/MVGC_result/coupling_resample_detrend_ERP(all_beta)';
    results_folder=fullfile(results_folder,subject);
    if ~exist(results_folder, 'dir')
        mkdir(results_folder);
    end
 save(fullfile(results_folder, 'TF_coupling.mat'), 'timedomain_G','-v7.3');   