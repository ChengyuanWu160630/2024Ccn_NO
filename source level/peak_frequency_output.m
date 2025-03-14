function[peak_freq]=peak_frequency_output(EEG,variable,foi,toi,A,B,C)
ft_data = eeglab2fieldtrip(EEG, 'preprocessing', 'none');
cfg = [];
switch variable
    case 'alpha'
    old_values=999999;
    new_values=0;       
        if toi(end)==1.3
            cfg.channel='E98';
        elseif toi(end)==1.35
            cfg.channel='E34';
        else 
            cfg.channel='E111';
        end
    case 'f'
    old_values=-999999;
    new_values=0;        
        if toi(end)>=1.24
            cfg.channel='E61';
        else
            cfg.channel='E111';
        end 
    case 'pwe'
        if toi(end)>0.5
        old_values=-999999;
        end
        
        if toi(end)==0.5
        old_values=999999;
        end    
        
        new_values=0;
        if toi(end)==0.5
            cfg.channel='E105';
        else
            cfg.channel='E39';
        end     
end
cfg.method = 'tfr';
cfg.output = 'pow';
cfg.foi = foi;
cfg.toi = [-0.6:0.01:-0.4,toi];
if max(foi)>13
    cfg.width = 7;
else
    cfg.width = 3;
end
cfg.keeptrials  = 'yes';
cfg.pad         =  'nextpow2';
subj_data_tf = ft_freqanalysis(cfg,ft_data);
cfg_baseline = [];
cfg_baseline.baseline = [-0.6 -0.4];
cfg_baseline.baselinetype = 'absolute';
subj_data_tf_blc = ft_freqbaseline(cfg_baseline, subj_data_tf);
slope=zeros(length(foi),1);
for i=1:length(foi)    
cfg_sel = [];
cfg_sel.latency = [toi];  
cfg_alpha = cfg_sel;
cfg_alpha.frequency = [foi(i)];
alpha_data = ft_selectdata(cfg_alpha, subj_data_tf_blc);
ALPHA=zscore(alpha_data.powspctrm);
        y=ALPHA;
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
                coefficients = mdl.Coefficients;
                new_values=coefficients.tStat(2); 
                slope(i)=coefficients.Estimate(2);
                if old_values>new_values
                    peak_freq=foi(i);
                    old_values=new_values;
                end               
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
                coefficients = mdl.Coefficients;
                new_values=coefficients.tStat(2);   
                slope(i)=coefficients.Estimate(2);               
                if old_values<new_values
                    peak_freq=foi(i);
                    old_values=new_values;
                end                
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
                coefficients = mdl.Coefficients;
                new_values=coefficients.tStat(2);
                slope(i)=coefficients.Estimate(2);
             if toi(end)>0.5   
                if old_values<new_values
                    peak_freq=foi(i);
                    old_values=new_values;
                end
             else
                if old_values>new_values
                    peak_freq=foi(i);
                    old_values=new_values;
                end
             end
        end
end
