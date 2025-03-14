function LME_calculation(theta,A,B,C,results_folder)
M = size(theta{1}.data, 1); % Number of electrodes
T = size(theta{1}.data, 2); % Number of time points
S = length(theta); % Number of subjects
t_values = zeros(M, T);
p_values = zeros(M, T);
masks = zeros(M, T);
y = cell(S, 1);
alpha = [];
f = [];
congruency = [];
pwe = [];
alpha_f = [];
alpha_congruency = [];
alpha_pwe = [];
subject_ids = [];
name3var={'alpha','f','pwe'};
receptor3var={A,B,C};
for testvar=1:length(name3var)
    variable=name3var{testvar};
    VAR=receptor3var{testvar};
    alpha = [];
    f = [];
    congruency = [];
    pwe = [];
    alpha_f = [];
    alpha_congruency = [];
    alpha_pwe = [];
    subject_ids = [];
for k = 1:S
    N_trials = size(theta{k}.data, 3); 
    alpha = [alpha; double(VAR{k}(1, 1:N_trials)')];
    f = [f; double(VAR{k}(2, 1:N_trials)')];
    congruency = [congruency; double(VAR{k}(3, 1:N_trials)')];
    pwe = [pwe; double(VAR{k}(4, 1:N_trials)')];
    alpha_f = [alpha_f; double(VAR{k}(5, 1:N_trials)')];
    alpha_congruency = [alpha_congruency; double(VAR{k}(6, 1:N_trials)')];
    alpha_pwe = [alpha_pwe; double(VAR{k}(7, 1:N_trials)')];
    subject_ids = [subject_ids; k * ones(N_trials, 1)];
end
if ~exist(fullfile(results_folder, variable), 'dir')
    mkdir(fullfile(results_folder, variable));
end
for i = 1:M
    for j = 1:T
        y = [];
        for k = 1:S
            N_trials = size(theta{k}.data, 3); 
            y = [y; squeeze(double(theta{k}.data(i, j, 1:N_trials)))];
        end
   switch variable
    case 'alpha'
        tbl = table(y, alpha, f, pwe, congruency, alpha_f, alpha_congruency, alpha_pwe, subject_ids, ...
                    'VariableNames', {'EEG', 'alpha', 'f', 'pwe', 'congruency', 'alpha_f', 'alpha_congruency', 'alpha_pwe', 'Subject'});
        simple_formula = 'EEG ~ 1 + alpha + f + pwe + congruency + alpha_f + alpha_congruency + alpha_pwe + (1|Subject)';
        complex_formula = 'EEG ~ 1 + alpha + f + pwe + congruency + alpha_f + alpha_congruency + alpha_pwe + (1 + alpha | Subject)';
    case 'f'
        tbl = table(y, f, alpha, pwe, congruency, alpha_f, alpha_congruency, alpha_pwe, subject_ids, ...
                    'VariableNames', {'EEG', 'f', 'alpha', 'pwe', 'congruency', 'alpha_f', 'alpha_congruency', 'alpha_pwe', 'Subject'});
        simple_formula = 'EEG ~ 1 + f + alpha + pwe + congruency + alpha_f + alpha_congruency + alpha_pwe + (1|Subject)';
        complex_formula = 'EEG ~ 1 + f + alpha + pwe + congruency + alpha_f + alpha_congruency + alpha_pwe + (1 + f | Subject)';
    case 'pwe'
        tbl = table(y, pwe, alpha, f, congruency, alpha_f, alpha_congruency, alpha_pwe, subject_ids, ...
                    'VariableNames', {'EEG', 'pwe', 'alpha', 'f', 'congruency', 'alpha_f', 'alpha_congruency', 'alpha_pwe', 'Subject'});
        simple_formula = 'EEG ~ 1 + pwe + alpha + f + congruency + alpha_f + alpha_congruency + alpha_pwe + (1|Subject)';
        complex_formula = 'EEG ~ 1 + pwe + alpha + f + congruency + alpha_f + alpha_congruency + alpha_pwe + (1 + pwe | Subject)';
    end
simple_mdl = fitlme(tbl, simple_formula);
complex_mdl = fitlme(tbl, complex_formula);
comparsion=compare(simple_mdl, complex_mdl);
pValue=comparsion.pValue(2);
ir_chi=comparsion.LRStat(2);
if pValue < 0.05
    selected_mdl = complex_mdl;
else
    selected_mdl = simple_mdl;
end
coefficients = selected_mdl.Coefficients;
t_values(i, j) = coefficients.tStat(2); 
p_values(i, j) = coefficients.pValue(2); 
masks(i, j) = p_values(i, j) < 0.05;
chi_value(i,j)=ir_chi;
choice(i,j)=pValue < 0.05;
    end
end
save(fullfile(results_folder, variable, 'results.mat'), 't_values', 'p_values', 'masks','choice','chi_value');
end