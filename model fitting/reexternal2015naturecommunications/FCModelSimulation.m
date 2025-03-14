%The flexible control model, originally implemented by Jiefeng Jiang 
%(jiefeng-jiang@uiowa.edu). Feel free to use/modify it for educational and/or 
%research purposes (at your own risk ;)). If you use this model in your 
%work, please cite our papers:
%
% 1. Jiang, J., Beck, J., Heller, K., & Egner, T. (2015). An insula-
% frontostriatal network mediates flexible cognitive control by adaptively 
% predicting changing control demands. Nature Communications, 
% DOI: 10.1038/ncomms9165 
%
% 2. Jiang, J., Heller, K., Egner, T.(2014). Bayesian modeling of flexible 
% cognitive control. Neuroscience and Biobehavioral Reviews, 46: 30-43.
%
%Input parameters:
%c: a sequence of (in)congruency (1 for incongruent, 0 for congruent, 
%   -1 for trials you'd like to skip).
%RTs: (transfomration of) behavioral data with normal distributions. If set
%     empty, this model only learns from c. 
%paras: parameters that maps predicted conflict level to RT. The script
%     will ignore this input if RTs is empty. 
%
%Output parameters:
%resultModel: an output model that records most of the information in
%            simulation.
%vEstimate: a vector that stores trial-by-trial estimates of learning rate. 
%
%pEstimate: a vector that stores trial-by-trial predicted conflict model.
%kEstimate: a vector that stores trial-by-trial k estimates.
%stdEstimate: a vector that stores trial-by-trial uncertainty of estimation.

function [resultModel, vEstimate, pEstimate, kEstimate, stdEstimate] = FCModelSimulation(c, RTs, paras)

model = FCModelInitialParameters;

resultModel = model;
pValues = model.para.pValue;
vValues = 2 .^(model.para.vValue);
mat = model.posterior;
tmp = mat;
dim = size(mat);
vEstimate = zeros(size(c));
pEstimate = zeros(size(c));
kEstimate = zeros(size(c));
stdEstimate = zeros(size(c));
jointPV = cell(1, length(c));

for trialCount = 2 : length(c)
    congruency = c(trialCount);
    lastCongruency = c(trialCount - 1);
    if mod(trialCount, 100) == 0  
        fprintf('Simulating trial %d\n', trialCount);
    end
    
    %Start making predictions of latent variables
    transVK = model.transVK;
    for j = 1 : dim(2)
        for k = 1 : dim(3)
            y = transVK(j, :, k)';
            x = mat(:, :, k);
            tmp(:, j, k) = x * y;
        end
    end
    mat = tmp;
    mat = mat / sum(sum(sum(mat)));
    
    %The "propagation" step
    for i = 1 : dim(1)
        for j = 1 : dim(2)
            x = model.transPV(i, :, j);
            y = reshape(mat(:, j, :), [dim(1) dim(3)]);
            tmp(i, j, :) = x * y;
            
        end
    end
    mat = tmp;
    
    %normalization
    mat = mat / sum(sum(sum(mat)));
    
    tmp(:, :, :) = 0;
    %updating the posterior distribution(fi+1->fi+0.5的转换，通过fi+0.5，vi+1，k的后验联合分布的插值（mat），比如说，step2.2有0.2的step3，0.8的step2)
    for i = 1 : dim(1)
        tPValues = pValues(i);
        for j = 1 : dim(2)
            weight = 1 / vValues(j);
            prediction = tPValues * (1 - weight) +  lastCongruency * weight;
            targetP = (prediction - model.para.pMin) / model.para.pStep + 1;
            if floor(targetP) < dim(1)
                weight1 = floor(targetP + 1) - targetP;
                weight2 = 1 - weight1;
                tmp(floor(targetP), j, :) = tmp(floor(targetP), j, :) + mat(i, j, :) * weight1;
                tmp(floor(targetP) + 1, j, :) = tmp(floor(targetP) + 1, j, :) + mat(i, j, :) * weight2;
            else
                tmp(floor(targetP), j, :) = tmp(floor(targetP), j, :) + mat(i, j, :);
            end
        end
    end
    mat = tmp;
    
    %normalization
    mat = mat / sum(sum(sum(mat)));
    
    %get estimates
    est = sum(mat, 3);
    ev = sum(est, 1);
    ep = squeeze(sum(est, 2));
    
    %using mean as estimate
    vEstimate(trialCount) = sum((1 ./ vValues) .* ev .* model.para.vdV);
    estK = sum(sum(mat, 2), 1);
    estK = reshape(estK, size(model.para.kValue));
    kEstimate(trialCount) = sum(model.para.kValue .* estK);
    pEstimate(trialCount) = sum(model.para.pValue .* ep');
    estimate = pEstimate(trialCount);
     
    sSum = 0;
    %updating posterior distribution & calculating prediction uncertainty 
    %(estimated as the std of predicted conflict)
    for i = 1 : length(ep)
        sSum = sSum + ep(i) * (estimate - model.para.pValue(i)) ^ 2;
    end
    stdEstimate(trialCount) = sqrt(sSum);
    
    %if congruency is less than 0, it means this trial is not used for
    %simulation.
    if congruency < 0
        %save joint distribution of probability and volatility as previous
        %trial
        vEstimate(trialCount) = vEstimate(trialCount - 1);
        pEstimate(trialCount) = pEstimate(trialCount - 1);
        if isSaveResults
            jointPV{trialCount} = est;
        end
        continue;
    end
    
    %updating belief when congruency and RT (if applicable) are observed.
    if RTs(trialCount) < 0
        for i = 1 : dim(1)
            tPValues = pValues(i);
            for j = 1 : dim(2)
                prediction = tPValues;
                factor = 1 - abs(congruency - prediction);
                mat(i, j, :) = mat(i, j, :) * factor;
            end
        end
    else
        for i = 1 : dim(1)
            tPValues = pValues(i);
            for j = 1 : dim(2)
                prediction = tPValues;
                if congruency > 0
                    rtMeanEst = paras(4) * prediction + paras(5);
                    factor1 = exp(-(RTs(trialCount) - rtMeanEst) * (RTs(trialCount) - rtMeanEst)) / 2 / paras(6) / paras(6);
                else
                    rtMeanEst = paras(1) * prediction + paras(2);
                    factor1 = exp(-(RTs(trialCount) - rtMeanEst) * (RTs(trialCount) - rtMeanEst)) / 2 / paras(3) / paras(3);
                end
                factor = (1 - abs(congruency - prediction)) * factor1;
                mat(i, j, :) = mat(i, j, :) * factor;
            end
        end
    end
    
    %normalization
    mat = mat / sum(sum(sum(mat)));
    
end
%set the prediction of the first trial to 0.5, as no privous info available
pEstimate(1) = 0.5;

%copy results
resultModel.posterior = mat;
