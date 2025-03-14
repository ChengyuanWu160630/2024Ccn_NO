function[sourceDiff]=wavelet_source_dpss_detect(subject,variable,toi,foi)
datapath=['/home/ChenQi_01/Downloads/ningbo_data/ningbo_behavior/behavior_result/',subject];
load(fullfile(datapath,[subject,'.mat']));
del_list=[];
for i=1:length(RTs)
    del=find(RTs{i}==-1)+80*(i-1);
    del_list=[del_list,del];
end
congruency=[c{1:8}];
for i=1:8
    load(strcat(datapath,'\sim',num2str(i),'.mat'))
    alpha((80*i-79):(80*i))=vEstimate;
    f((80*i-79):(80*i))=pEstimate;
    uncertainty((80*i-79):(80*i))=stdEstimate;
end
fn=[subject,'.set'];
data_path=['/home/Downloads/ningbo_data/finished'];
EEG = pop_loadset('filename',fn,'filepath',data_path);
firstdel=EEG.manual_refuse;
seconddel=EEG.auto_refuse;
EEGdel={firstdel,seconddel};
%% remove epoch
EEG_R=EEG;
f_R=f;
alpha_R=alpha;
congruency_R=congruency;
uncertainty_R=uncertainty;
[~,~,~,f,alpha,pwe,congruency,EEG]=EEG_align_behavior(f,alpha,congruency,uncertainty,EEG,EEGdel,del_list,0);
[~,~,~,R_f,R_alpha,R_pwe,~]=EEG_align_behavior(f_R,alpha_R,congruency_R,uncertainty_R,EEG_R,EEGdel,del_list,1);
EEG = pop_epoch( EEG, {  '2   '  }, [-1           3.002], 'newname', 'EEProbe continuous data epochs', 'epochinfo', 'yes');
A=[R_alpha;f;congruency;pwe;alpha.*f;alpha.*congruency;alpha.*pwe];
B=[alpha;R_f;congruency;pwe;alpha.*f;alpha.*congruency;alpha.*pwe];
C=[alpha;f;congruency;R_pwe;alpha.*f;alpha.*congruency;alpha.*pwe];
[peak_freq]=peak_frequency_output(EEG,variable,foi,toi,A,B,C);
%% leadfield
templateheadmodel = '/home/Downloads/egi/vol.mat';
load(templateheadmodel); 
elecPath='/home/Downloads/egi/elec_aligned.mat';
load(elecPath);
cfg                 = [];
cfg.elec =  elec_aligned;    
cfg.headmodel = vol;        
cfg.reducerank      = 3;
cfg.channel         = 'all';
cfg.grid.resolution = 1;   
cfg.grid.unit       = 'cm';
cfg.normalize       = 'yes';
leadfield = ft_prepare_leadfield(cfg);
%% data preparation
ft_data = eeglab2fieldtrip(EEG, 'preprocessing', 'none');
cfg = [];                                           
cfg.toilim = [-0.75 -0.25];         
dataPre = ft_redefinetrial(cfg, ft_data);
cfg.toilim = [toi-0.25 toi+0.25];                       
dataPost = ft_redefinetrial(cfg, ft_data);
design = [ones(1,length(dataPre.trial)) ones(1,length(dataPost.trial))*2]; 
%% frequency analysis 
dataAll = ft_appenddata([], dataPre, dataPost); 
cfg=[];
cfg.method      = 'mtmfft';
cfg.output      = 'powandcsd';  
cfg.foilim      = [peak_freq peak_freq];     
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 3;
cfg.keeptrials  = 'yes';     
cfg.pad         =  'nextpow2';
freqALL = ft_freqanalysis(cfg, dataAll);
%% common spatial filter 
cfg=[];
cfg.method      = 'dics';
cfg.sourcemodel = leadfield;         
cfg.headmodel   = vol;         
cfg.frequency   = freqALL.freq;
cfg.dics.keepfilter  = 'yes';    
cfg.dics.realfilter   = 'yes';
cfg.dics.lambda       = '10%';
source = ft_sourceanalysis(cfg, freqALL);
% project all trials through common spatial filter 
cfg=[];
cfg.method      = 'dics';
cfg.sourcemodel = leadfield;       
cfg.headmodel   = vol;        
cfg.sourcemodel.filter = source.avg.filter; 
cfg.frequency   = freqALL.freq;
cfg.rawtrial    = 'yes';      
source = ft_sourceanalysis(cfg, freqALL); 
% extract single trial result for source power 
AA = find(design==1); % find trial numbers belonging to baseline
BB = find(design==2); % find trial numbers belonging to post stimulus
sourcePre = source;
sourcePre.trial(BB) = [];
sourcePre.cumtapcnt(BB) = [];
sourcePre.df = length(AA);
cfg.keeptrials       = 'yes';
sourcePre = ft_sourcedescriptives(cfg, sourcePre); 
sourcePost=source;
sourcePost.trial(AA) = [];
sourcePost.cumtapcnt(AA) = [];
sourcePost.df = length(BB);
sourcePost = ft_sourcedescriptives(cfg, sourcePost);
sourceDiff = ft_sourcedescriptives([], sourcePost); 
%% fitting GLM
voxel_activity=zeros(length(sourceDiff.avg.pow),EEG.trials);
voxel_activity2=zeros(length(sourceDiff.avg.pow),EEG.trials);
for i=1:EEG.trials
sourceDiff.avg.pow = (sourcePost.trial(i).pow - sourcePre.trial(i).pow) ./ sourcePre.trial(i).pow;
sourceDiff.avg.pow2 = (sourcePost.trial(i).pow - sourcePre.trial(i).pow) ;
voxel_activity(:,i)=sourceDiff.avg.pow;
voxel_activity2(:,i)=sourceDiff.avg.pow2;
    if i==1
        nanIndicesPre = ~isnan(sourcePost.trial(i).pow);
    end
end
GLMvoxel=voxel_activity(nanIndicesPre,:);
GLMvoxel2=voxel_activity2(nanIndicesPre,:);

slope_values=zeros(size(GLMvoxel,1),1);
slope_tvalue=zeros(size(GLMvoxel,1),1);
for k=1:size(GLMvoxel,1)
        y=zscore(GLMvoxel(k,:)');
        switch variable
            case 'alpha'               
                alpha = A(1,:)';
                f = A(2,:)';
                congruency = A(3,:)';
                pwe = A(4,:)';
                alpha_f = A(5,:)';
                alpha_congruency = A(6,:)';
                alpha_pwe = A(7,:)';
                tbl = table(y, alpha,f,pwe,congruency,alpha_f, alpha_congruency, alpha_pwe, ...
                'VariableNames', {'EEG', 'alpha', 'f','pwe','congruency','alpha_f', 'alpha_congruency', 'alpha_pwe'});      
                formula = 'EEG ~ alpha +f+pwe +congruency+ alpha_f+ alpha_congruency + alpha_pwe';       
                mdl = fitlm(tbl, formula);
            case 'f'
                alpha = B(1,:)';
                f = B(2,:)';
                congruency = B(3,:)';
                pwe = B(4,:)';
                alpha_f = B(5,:)';
                alpha_congruency = B(6,:)';
                alpha_pwe = B(7,:)';
                tbl = table(y, f, alpha,pwe,congruency,alpha_f, alpha_congruency, alpha_pwe, ...
                'VariableNames', {'EEG', 'f', 'alpha','pwe','congruency','alpha_f', 'alpha_congruency', 'alpha_pwe'});      
                formula = 'EEG ~ f+alpha+pwe + congruency+alpha_f + alpha_congruency + alpha_pwe';       
                mdl = fitlm(tbl, formula);
            case 'pwe'
                alpha = C(1,:)';
                f = C(2,:)';
                congruency = C(3,:)';
                pwe = C(4,:)';
                alpha_f = C(5,:)';
                alpha_congruency = C(6,:)';
                alpha_pwe = C(7,:)';
                tbl = table(y, pwe, alpha,f,congruency,alpha_f, alpha_congruency, alpha_pwe, ...
                'VariableNames', {'EEG', 'pwe','alpha','f','congruency', 'alpha_f', 'alpha_congruency', 'alpha_pwe'});      
                formula = 'EEG ~ pwe +alpha+f +congruency+ alpha_f + alpha_congruency + alpha_pwe';       
                mdl = fitlm(tbl, formula);
        end
        coefficients = mdl.Coefficients;
        slope_values(k)=coefficients.Estimate(2);
        slope_tvalue(k)=coefficients.tStat(2);
end
sourceDiff.avg.pow(nanIndicesPre)=slope_tvalue;
slope_values=zeros(size(GLMvoxel2,1),1);
slope_tvalue=zeros(size(GLMvoxel,1),1);
for k=1:size(GLMvoxel2,1)
        y=zscore(GLMvoxel2(k,:)');
        switch variable
            case 'alpha'               
                alpha = A(1,:)';
                f = A(2,:)';
                congruency = A(3,:)';
                pwe = A(4,:)';
                alpha_f = A(5,:)';
                alpha_congruency = A(6,:)';
                alpha_pwe = A(7,:)';
                tbl = table(y, alpha,f,pwe,congruency,alpha_f, alpha_congruency, alpha_pwe, ...
                'VariableNames', {'EEG', 'alpha', 'f','pwe','congruency','alpha_f', 'alpha_congruency', 'alpha_pwe'});      
                formula = 'EEG ~ alpha +f+pwe +congruency+ alpha_f+ alpha_congruency + alpha_pwe';       
                mdl = fitlm(tbl, formula);
            case 'f'
                alpha = B(1,:)';
                f = B(2,:)';
                congruency = B(3,:)';
                pwe = B(4,:)';
                alpha_f = B(5,:)';
                alpha_congruency = B(6,:)';
                alpha_pwe = B(7,:)';
                tbl = table(y, f, alpha,pwe,congruency,alpha_f, alpha_congruency, alpha_pwe, ...
                'VariableNames', {'EEG', 'f', 'alpha','pwe','congruency','alpha_f', 'alpha_congruency', 'alpha_pwe'});      
                formula = 'EEG ~ f+alpha+pwe + congruency+alpha_f + alpha_congruency + alpha_pwe';       
                mdl = fitlm(tbl, formula);
            case 'pwe'
                alpha = C(1,:)';
                f = C(2,:)';
                congruency = C(3,:)';
                pwe = C(4,:)';
                alpha_f = C(5,:)';
                alpha_congruency = C(6,:)';
                alpha_pwe = C(7,:)';
                tbl = table(y, pwe, alpha,f,congruency,alpha_f, alpha_congruency, alpha_pwe, ...
                'VariableNames', {'EEG', 'pwe','alpha','f','congruency', 'alpha_f', 'alpha_congruency', 'alpha_pwe'});      
                formula = 'EEG ~ pwe +alpha+f +congruency+ alpha_f + alpha_congruency + alpha_pwe';       
                mdl = fitlm(tbl, formula);
        end
        coefficients = mdl.Coefficients;
        slope_values(k)=coefficients.Estimate(2);
        slope_tvalue(k)=coefficients.tStat(2);
end
sourceDiff.avg.pow2(nanIndicesPre)=slope_tvalue;