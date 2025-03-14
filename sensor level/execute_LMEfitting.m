cd('F:\eeg_ningbo\eeg_data\finished')
files=dir('*.set');
document={files.name};
theta = cell(length(document), 1);
ALPHA = cell(length(document), 1);
beta1 = cell(length(document), 1);
beta2 = cell(length(document), 1);
EEG_path='F:\eeg_ningbo\eeg_data\finished';
behavior_path='F:\eeg_ningbo\ningbo_behavior\behavior_result';
A = cell(length(document), 1);
B= cell(length(document), 1);
C = cell(length(document), 1);
D = cell(length(document), 1);
E = cell(length(document), 1);
numThreads = 2;
if isempty(gcp('nocreate'))
     parpool(numThreads, 'IdleTimeout', Inf);
end
parfor i=1:length(document)
    subject=document{i}(1:3);
    [theta{i},ALPHA{i},beta1{i},beta2{i},A{i},B{i},C{i},D{i},E{i}] = TF_LME_sti(subject,behavior_path,EEG_path); % A has the residuals for alpha
end
delete(gcp('nocreate'));
numThreads = 2; 
results_folder1='F:\results\LME\theta';
results_folder2='F:\results\LME\ALPHA';
if isempty(gcp('nocreate'))
    parpool(numThreads, 'IdleTimeout', Inf);
end
f1 = parfeval(@LME_calculation, 0, theta, A, B, C,  results_folder1);
f2 = parfeval(@LME_calculation, 0, ALPHA, A, B, C,  results_folder2);
wait(f1);
wait(f2);
delete(gcp('nocreate'));