function alpha_MVGC_granger(subject)
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
group5 ={'E112','E106','E105','E80','E79','E6','E37','E16'};%E111 E87 delete
all_electrodes = [group1, group2, group3, group4,group5];
EEG = pop_select( EEG, 'channel',all_electrodes);
EEG = pop_resample( EEG, 200);
groups = {group1, group2, group3, group4};
pairs = {};
for i = 1:length(groups)
    for j = 1:length(groups)
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
% substractERP
mean_EEG=squeeze(mean(EEG.data,3));
for i=1:size(EEG.data,3)
    EEG.data(:,:,i)=EEG.data(:,:,i)-mean_EEG;

end
%% Parameters
regmode   = 'LWR';   
fs        = 200;   
fres      = 1000;    
morder    = 15;    
win= 100;
tstat     = 'F';   
%% Vertical regression GC calculation
glanger_time=round([-0.6:0.01:2.7]*1000);
frequency_G=zeros([2,fres+1,length(glanger_time),size(pairs,1)],'single');
for elec=1:size(pairs,1)
    idx=index_pairs(elec,:);
    X=EEG.data(idx,:,:);
    for e = 1:length(glanger_time)
        j = find(EEG.times==glanger_time(e));
        slidewin=X(:,j-win/2:j+win/2,:);
        for T=1:size(slidewin,3)
            slidewin(1,:,T)=zscore(detrend(squeeze(slidewin(1,:,T))));
            slidewin(2,:,T)=zscore(detrend(squeeze(slidewin(2,:,T))));
        end
        try
            [A,SIG] = tsdata_to_var(slidewin,morder,regmode);
        catch
            continue
        end
        info = var_info(A,SIG);
        if isempty(fres)
        fres = 2^nextpow2(info.acdec); 
        fprintf('\nfrequency resolution auto-calculated as %d (increments ~ %.2gHz)\n',fres,fs/2/fres);
        end        
        if fres > 20000 
        fprintf(2,'\nWARNING: large frequency resolution = %d - may cause computation time/memory usage problems\nAre you sure you wish to continue [y/n]? ',fres);
        istr = input(' ','s'); if isempty(istr) || ~strcmpi(istr,'y'); fprintf(2,'Aborting...\n'); return; end
        end
        ptic('\n*** var_to_spwcgc... ');
        freq_G = var_to_spwcgc(A,SIG,fres);
        ptoc;
        frequency_G(1,:,e,elec)=squeeze(freq_G(2,1,:));
        frequency_G(2,:,e,elec)=squeeze(freq_G(1,2,:));
    end
end
results_folder='/home/ChenQi_01/Downloads/method_wavelet/direction_connection/granger/mvgc/MVGC_result/resample200(step10ms)_erp_detrend';
    results_folder=fullfile(results_folder,subject);
    if ~exist(results_folder, 'dir')
        mkdir(results_folder);
    end
 save(fullfile(results_folder, 'TF_granger.mat'), 'frequency_G','-v7.3');   