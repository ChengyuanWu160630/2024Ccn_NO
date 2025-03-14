%Initialize model parameters
function BDModel = FCModelInitialParameters

%k: distribution factor of volatility, logged
%v: volatility, logged
%p: estimate of the probability that the incoming trial is incongruent
BVpara.pMin = 0;
BVpara.pMax = 1;
BVpara.pStep = .02;
BVpara.pValue = BVpara.pMin: BVpara.pStep: BVpara.pMax;

%change the range of volatility/learning rate based on your task design
BVpara.vMin = 1;
BVpara.vMax = 160; 
BVpara.vStep = (BVpara.vMax - BVpara.vMin) / 79; 
BVpara.vValue = BVpara.vMin : BVpara.vStep : BVpara.vMax;
BVpara.vdV = zeros(size(BVpara.vValue));
BVpara.vdV(1) = (BVpara.vValue(2) - BVpara.vValue(1)) / 2;
BVpara.vdV(end) = (BVpara.vValue(end) - BVpara.vValue(end - 1)) / 2;
for i = 2 : length(BVpara.vValue) - 1
    BVpara.vdV(i) = (BVpara.vValue(i + 1) - BVpara.vValue(i - 1)) / 2;
end
BVpara.vdV = BVpara.vdV / BVpara.vdV(2);
BVpara.vValue = log(BVpara.vValue) / log(2);

%gives you a full range of meta-volatility
BVpara.kMin = 0.04;
BVpara.kMax = 0.96;
BVpara.kStep = 0.04;
BVpara.kValue = BVpara.kMin : BVpara.kStep : BVpara.kMax;

BDModel.para = BVpara;

%posterior distribution, no prior knowledge, which means posteior
%distribution is initialized as a uniform distribution. The matrix is
%organized in the order of p_{i+1}, v_{i+1}, k_{i+1}
posterior = ones([length(BVpara.pValue) length(BVpara.vValue) length(BVpara.kValue)]);
posterior = posterior / (length(BVpara.pValue) * length(BVpara.vValue) * length(BVpara.kValue));
BDModel.posterior = posterior;

%initialize translation distributions

%initialize transVK p(v_{i+1}|v_i,k_i) as a Gaussian distribution
%the dimensions of this 3D matrix are organized in the
%order of v_{i+1}, v_i, k_i
transVK = zeros(length(BVpara.vValue), length(BVpara.vValue), length(BVpara.kValue));
for k = 1 : length(BVpara.kValue)
    for i = 1 : length(BVpara.vValue)
        for j = 1 : length(BVpara.vValue)
           transVK(i, j, k) = (1 - BVpara.kValue(k)) / length(BVpara.vValue);
        end
        transVK(i, i, k) = transVK(i, i, k) + BVpara.kValue(k);
    end
end

for j = 1 : length(BVpara.vValue)
    for k = 1 : length(BVpara.kValue)
        col = transVK(:, j, k);
        sumCol = sum(col);
        if sumCol > 0
            col = col / sumCol;
        else
            col = ones(length(BVpara.vValue), 1) / length(BVpara.vValue);
        end
        transVK(:, j, k) = col(:, :);
    end
end


BDModel.transVK = transVK;

%initialize transPV transPV p(p_{i+1}|p_i,v_{i+1}) as a Beta distribution
%the dimensions of this 3D matrix are organized in the
%order of p_{i+1}, p_i, v_{i+1}
%the logrithm computations reduce magnitudes of intermediate results and in turn the chance of having Inf or NaN results 
transPV = zeros(length(BVpara.pValue), length(BVpara.pValue), length(BVpara.vValue));
for k = 1 : length(BVpara.vValue)
    %v_ip1 = exp(BVpara.vValue(k));
    v_ip1 = 2 ^ BVpara.vValue(k) + 2;
    gammaSum = gammaln(v_ip1);
    for j = 1 : length(BVpara.pValue)
        p_i = BVpara.pValue(j);
        %if the sharpness (alpha + beta) is no greater than 2 (meaning mode is not available), center the mean to p_i, otherwise center the mode to p_i; center the sharpness as v_{i+1}
        if v_ip1 <= 2
            alpha = p_i * v_ip1;
        else
            alpha = p_i * (v_ip1 - 2) + 1;
        end
        beta = v_ip1 - alpha;
        x = gammaln(alpha) + gammaln(beta) - gammaSum;
        for i = 1 : length(BVpara.pValue)
            p_ip1 = BVpara.pValue(i);
            if (p_ip1 == 0 || p_ip1 == 1)
                transPV(i, j, k) = 0;
            else
                transPV(i, j, k) = exp(x + log(p_ip1) * (alpha - 1) + log(1 - p_ip1) * (beta - 1));
            end
        end
    end
end

for j = 1 : length(BVpara.pValue)
    for k = 1 : length(BVpara.vValue)
        col = transPV(:, j, k);
        sumCol = sum(col);
        if sumCol > 0
            col = col / sum(col);
        else
            col = ones(length(BVpara.pValue), 1) / length(BVpara.pValue);
        end
        transPV(:, j, k) = col(:);
    end
end


BDModel.transPV = transPV;    
    
    



