cd('F:\eeg_ningbo\ningbo_behavior\behavior_result');
files=dir();
for d=1:length(files)
path=['F:\eeg_ningbo\ningbo_behavior\behavior_result\',files(d).name];
loadfile=[files(d).name,'.mat'];
cd(path);
load(loadfile);
    for i=1:8
    RT_1=(RTs{i}>0);
    RTs{i}(RT_1)=RTs{i}(RT_1)*1000;
    FCModelRT(RTs{i}, c{i},i)
    end
end
