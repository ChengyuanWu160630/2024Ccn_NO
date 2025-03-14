function MVGC_coupling_theta_alpha(subject)
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

% define significant electrodes
group1 = {'E32','E26','E28','E29','E33','E34','E35','E36','E39','E40','E45'};
group2 = {'E86','E92','E97','E98','E93','E87'};
group3 = {'E1','E2','E8','E9','E122','E123','E124','E117','E116','E110','E111','E109','E108','E104','E103'};
group4 = {'E47','E42','E52','E53','E60','E66','E54','E61','E67','E62','E72','E78','E77','E84','E85','E90','E91'};
group5 ={'E112','E106','E105','E80','E79','E6','E37','E16'};
all_electrodes = [group1, group2, group3, group4,group5];
EEG = pop_select( EEG, 'channel',all_electrodes);

groups = {group1, group2, group3, group4, group5};
pairs = {};
for i = 1:length(groups)
    j=5;
        if j > i
            g1 = groups{i};
            g2 = groups{j};          
            for k = 1:length(g1)
                for l = 1:length(g2)
                    pairs{end+1, 1} = g1{k};
                    pairs{end, 2} = g2{l};
                end
            end
        end  
end
pairs = reshape(pairs, [], 2);
labelslist = {EEG.chanlocs.labels};
index_pairs = zeros(size(pairs));
for i = 1:size(pairs, 1)
    for j = 1:size(pairs, 2)
        index_pairs(i, j) = find(strcmp(labelslist, pairs{i, j}));
    end
end
firstdel=EEG.manual_refuse;
seconddel=EEG.auto_refuse;
EEGdel={firstdel,seconddel};
%% remove epoch
[~,~,~,~,~,~,~,EEG]=EEG_align_behavior(f,alpha,congruency,uncertainty,EEG,EEGdel,del_list,0);
EEG = pop_epoch( EEG, {  '2   '  }, [-1           3.002], 'newname', 'EEProbe continuous data epochs', 'epochinfo', 'yes');
%% substitute data
ft_data = eeglab2fieldtrip(EEG, 'preprocessing', 'none');
load(fullfile('/home/ChenQi_01/Downloads/method_wavelet/direction_connection/granger/nonpara/coupling/FRE_nopading',subject,'FRE.mat'));
labels = ft_data.label;
indices = zeros(1, length(group5));
for i = 1:length(group5)
    idx = find(ismember(labels, group5{i}));    
    if ~isempty(idx)
        indices(i) = idx;
    else
        indices(i) = 1000;  
    end
end
ALPHA_SPM=squeeze(alpha_data.powspctrm);
theta_SPM=squeeze(theta_data.powspctrm);
ALPHA_SPM(isnan(ALPHA_SPM))=0;
theta_SPM(isnan(theta_SPM))=0;
for i=1:length(ft_data.trial)
ft_data.trial{i}=squeeze(ALPHA_SPM(i,:,:));
ft_data.trial{i}(indices,:)=squeeze(theta_SPM(i,indices,:));
end
cfg=[];
cfg.latency=[-0.63 2.63];
ft_data=ft_selectdata(cfg,ft_data);
cfg=[];
cfg.avgoverrpt ='yes';
ft_data_avg=ft_selectdata(cfg,ft_data);

EEG = pop_epoch( EEG, {  '2   '  }, [-0.63           2.632], 'newname', 'EEProbe continuous data epochs', 'epochinfo', 'yes');
for i=1:length(ft_data.trial)
    EEG.data(:,:,i)=ft_data.trial{i}-ft_data_avg.trial{1};
%    EEG.data(:,:,i)=ft_data.trial{i};
end
EEG = pop_resample( EEG, 200);
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
    X=EEG.data(idx,:,:);

for e = 1:length(glanger_time)
    j = find(EEG.times==glanger_time(e));
    slidewin=X(:,j-win/2:j+win/2,:);
    for T=1:size(slidewin,3)
        slidewin(1,:,T)=zscore(detrend(squeeze(slidewin(1,:,T))));
        slidewin(2,:,T)=zscore(detrend(squeeze(slidewin(2,:,T))));
%         slidewin(1,:,T)=zscore(squeeze(slidewin(1,:,T)));
%         slidewin(2,:,T)=zscore(squeeze(slidewin(2,:,T)));
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

results_folder='/home/ChenQi_01/Downloads/method_wavelet/direction_connection/granger/mvgc/MVGC_result/coupling(resample_date_Z_ERP_detrend1023)';
    results_folder=fullfile(results_folder,subject);
    if ~exist(results_folder, 'dir')
        mkdir(results_folder);
    end
 save(fullfile(results_folder, 'TF_coupling.mat'), 'timedomain_G','-v7.3');   