data_path = 'F:\eeg_ningbo\ningbo_behavior\behavior';
save_path='F:\eeg_ningbo\ningbo_behavior\behavior_result';
cd(data_path);
files=dir('*.mat');
behaviordata={files.name};
for i=1:length(files)
    cd(data_path);
    load(behaviordata{i});
    [run,RTs,c,dealdata] = behavior_preprocess(run640);
    cd(save_path);
    if ~exist(fullfile(save_path, behaviordata{i}(1:3)), 'dir')
    mkdir(fullfile(save_path, behaviordata{i}(1:3)));
    end
    cd(behaviordata{i}(1:3));
    savename=[behaviordata{i}(1:3),'.mat'];
    save(savename,'run','c','RTs','dealdata');
    clear('run','c','RTs','dealdata');
end
