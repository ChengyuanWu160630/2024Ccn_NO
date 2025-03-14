files=dir('F:\eeg_ningbo\wavelet_data\finished\*.set');
document={files.name};
group = cell(length(document), 1);
variable='alpha';
toi=1.35;
foi=[8:0.5:12];
numThreads = 2;
% % 如果当前没有活动的并行池，创建一个
% if isempty(gcp('nocreate'))
%      parpool(numThreads, 'IdleTimeout', Inf);
% end
for i=1:length(document)
    subject=document{i}(1:3);
            [group{i}]=wavelet_sourceEEG_detect(subject,variable,toi,foi);
end
% delete(gcp('nocreate'));

results_folder='F:\eeg_ningbo\source\TFmethod\dicsbytrials\wavelet\detect_singlepeakfrequency_and_source\results\averagetime\realfilter\10lambda';
if ~exist(fullfile(results_folder, variable), 'dir')
    mkdir(fullfile(results_folder, variable));
end
save (fullfile(results_folder,variable, 'source_alpha.mat'), 'group')