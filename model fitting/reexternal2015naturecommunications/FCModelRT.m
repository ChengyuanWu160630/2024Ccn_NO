
function FCModelRT(RTs, c,i) % 

    RT = 1000 ./ RTs; 
    RTIdxC = (c < 0.5) & (RTs > 0);
    RTIdxI = (c > 0.5) & (RTs > 0);
    
    RTC = RT(RTIdxC)';
    RTI = RT(RTIdxI)';
    
    %find initial value
    [~, ~, pEstimate, ~, ~] = FCModelSimulation(c, -ones(size(c)), []);
    
    %estimate parameters. For congruent trials, RT~N(pEstimate * bC(2) + bC(1), sigma1)
    %For incongruent trials, RT~N(pEstimate * bI(2) + bI(1), sigma2)
    
    pC = pEstimate(RTIdxC);
    bC = glmfit(pC', RTC);
    r1 = RTC - bC(2) * pC' - bC(1);
    sigma1 = sqrt(sum(r1 .^ 2) / length(r1));
   % pI = pEstimate(RTIdx(RTIdxI));
    pI = pEstimate(RTIdxI);  % MJ changed!
    bI = glmfit(pI', RTI);
    r1 = RTI - bI(2) * pI' - bI(1);
    sigma2 = sqrt(sum(r1 .^ 2) / length(r1));
    
    paras = [bC(2) bC(1) sigma1 bI(2) bI(1) sigma2]
    
    count = 0;
    dis = 1;
    
    %the EM algorithm the iteratively update parameter estimates and model
    %estimates
    while dis > 1e-3 && count < 100
        count = count + 1;
      
        
       [~, vEstimate, pEstimate, kEstimate, stdEstimate] = FCModelSimulation(c, RT, paras); 
         
    
        oldParas = paras;
        %pC = pEstimate(RTIdx(RTIdxC));
        pC = pEstimate(RTIdxC); % MJ changed!
        bC = glmfit(pC', RTC);
        r1 = RTC - bC(2) * pC' - bC(1);
        sigma1 = sqrt(sum(r1 .^ 2) / length(r1));
       % pI = pEstimate(RTIdx(RTIdxI));
        pI = pEstimate(RTIdxI); % MJ changed!
        bI = glmfit(pI', RTI);
        r1 = RTI - bI(2) * pI' - bI(1);
        sigma2 = sqrt(sum(r1 .^ 2) / length(r1));
        paras = [bC(2) bC(1) sigma1 bI(2) bI(1) sigma2]
%         paras(isnan(paras))=0;
        %calculate the difference between old and new estimates of
        %parameters
        dis = sqrt(sum((oldParas - paras) .^ 2))
    end
    NAME="sim"+i+".mat";
    %save results 
    pEstimate(isnan(pEstimate))=0;
    vEstimate(isnan(vEstimate))=0;
    save(NAME, 'vEstimate', 'pEstimate', 'stdEstimate', 'kEstimate', 'paras','count');
