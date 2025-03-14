files=dir('/home/ChenQi_01/Downloads/ningbo_data/finished/*.set');
document={files.name};
numThreads = 9;
if isempty(gcp('nocreate'))
     parpool(numThreads, 'IdleTimeout', Inf);
end
parfor i=1:length(files)
    subject=document{i}(1:3);
    MVGC_coupling_beta(subject);
end
delete(gcp('nocreate'));