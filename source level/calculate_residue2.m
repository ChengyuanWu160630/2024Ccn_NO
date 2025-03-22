function[residuals]=calculate_residue2(pwe,alpha,f,congruency,uncertainty,copyF_U,a_PE,copyalpha,need_to_changeZ)
pwe_noneed_deleted=find(need_to_changeZ==1);
input_alpha=[];
input_f=[];
input_pwe=[];
input_num=[];
input_congruency=[];
input_uncertainty=[];
input_a_PE=[];
input_copy_alpha=[];
% NOZ_f=[];
% NOZ_U=[];
input_F_U=[];
for i=1:length(pwe_noneed_deleted)
    extract_num=pwe_noneed_deleted(i);
    if isempty(find([80:80:640]==extract_num, 1)) && ~isempty(find(pwe_noneed_deleted==(extract_num+1), 1))
        input_pwe(end+1)=pwe(extract_num);
        input_congruency(end+1)=congruency(extract_num);
        input_f(end+1)=f(extract_num+1);
        input_alpha(end+1)=alpha(extract_num+1);
        input_copy_alpha=copyalpha(extract_num+1);
        input_uncertainty(end+1)=uncertainty(extract_num+1);
        input_a_PE(end+1)=a_PE(extract_num);        
%         NOZ_f(end+1)=copyF(extract_num+1);
%         NOZ_U(end+1)=copyU(extract_num+1);    
        input_F_U(end+1)=copyF_U(extract_num+1);
        input_num(end+1)=extract_num;
    end
end
   %% pwe
    % 构造设计矩阵并包括截距
    %X = [ones(size(input_alpha))', input_alpha', input_f', input_congruency',(input_alpha.*input_f)',(input_alpha.*input_congruency)',(input_alpha.*input_congruency.*input_f)'];
    X = [ones(size(input_alpha))', input_alpha', input_f', input_congruency',(input_alpha.*input_f)',(input_alpha.*input_congruency)',(input_alpha.*input_pwe)'];
   
    % 用pwe对alpha，f和congruency进行多元线性回归
    b = regress(input_pwe', X);
    % 计算并存储pwe的残差
    yPredicted = X * b;
    residue_pwe= input_pwe - yPredicted';
    pwe(input_num)=residue_pwe;
    residuals.pwe =pwe;
    

   %% f
      % 构造设计矩阵并包括截距
    %X = [ones(size(input_alpha))', input_alpha', input_pwe', input_congruency',(input_alpha.*input_f)',(input_alpha.*input_congruency)',(input_alpha.*input_congruency.*input_f)'];
    X = [ones(size(input_alpha))', input_alpha', input_pwe', input_congruency',(input_alpha.*input_f)',(input_alpha.*input_congruency)',(input_alpha.*input_pwe)'];

    % 用pwe对alpha，f和congruency进行多元线性回归
    b = regress(input_f', X);
    % 计算并存储pwe的残差
    yPredicted = X * b;
    residue_f= input_f - yPredicted';
    f(input_num+1)=residue_f;
    residuals.f =f;
    
   %% alpha
      % 构造设计矩阵并包括截距
    %X = [ones(size(input_alpha))', input_pwe', input_f', input_congruency',(input_alpha.*input_f)',(input_alpha.*input_congruency)',(input_alpha.*input_congruency.*input_f)'];
    X = [ones(size(input_alpha))', input_pwe', input_f', input_congruency',(input_alpha.*input_f)',(input_alpha.*input_congruency)',(input_alpha.*input_pwe)'];
  
    % 用pwe对alpha，f和congruency进行多元线性回归
    b = regress(input_alpha', X);
    % 计算并存储pwe的残差
    yPredicted = X * b;
    residue_alpha= input_alpha - yPredicted';    
    alpha(input_num+1)=residue_alpha;
    residuals.alpha=alpha;

    residuals.input_num=input_num;
    assert(length(residuals.alpha) == length(residuals.pwe), '维度不相等！');

   %% alpha_pwe
      % 构造设计矩阵并包括截距
    %X = [ones(size(input_alpha))', input_pwe', input_f', input_congruency',(input_alpha.*input_f)',(input_alpha.*input_congruency)',(input_alpha.*input_congruency.*input_f)'];
    X = [ones(size(input_alpha))', input_alpha',input_pwe', input_f', input_congruency',(input_alpha.*input_f)',(input_alpha.*input_congruency)'];
    alpha_pwe=1:length(alpha);
    % 用pwe对alpha，f和congruency进行多元线性回归
    b = regress((zscore(input_a_PE.*input_copy_alpha))', X);
    % 计算并存储pwe的残差
    yPredicted = X * b;
    residue_alpha_pwe= zscore(input_a_PE.*input_copy_alpha)- yPredicted';    
    alpha_pwe(input_num+1)=residue_alpha_pwe;
    residuals.alpha_pwe=alpha_pwe;
    
       %% uncertainty
      % 构造设计矩阵并包括截距
    %X = [ones(size(input_alpha))', input_pwe', input_f', input_congruency',(input_alpha.*input_f)',(input_alpha.*input_congruency)',(input_alpha.*input_congruency.*input_f)'];
    X = [ones(size(input_alpha))',input_alpha', input_pwe', input_f', input_congruency',(input_alpha.*input_f)',(input_alpha.*input_congruency)',(input_alpha.*input_pwe)'];
  
    % 用pwe对alpha，f和congruency进行多元线性回归
    b = regress(input_uncertainty', X);
    % 计算并存储pwe的残差
    yPredicted = X * b;
    residue_uncertainty= input_uncertainty - yPredicted';    
    uncertainty(input_num+1)=residue_uncertainty;
    residuals.uncertainty=uncertainty;
    
   %% uncertainty*f
      % 构造设计矩阵并包括截距
    %X = [ones(size(input_alpha))', input_pwe', input_f', input_congruency',(input_alpha.*input_f)',(input_alpha.*input_congruency)',(input_alpha.*input_congruency.*input_f)'];
    X = [ones(size(input_alpha))',input_alpha', input_pwe', input_f', input_congruency',(input_alpha.*input_f)',(input_alpha.*input_congruency)',(input_alpha.*input_pwe)',input_uncertainty'];
  
    % 用pwe对alpha，f和congruency进行多元线性回归
    b = regress(input_F_U', X);
    % 计算并存储pwe的残差
    yPredicted = X * b;
    residue_uncertainty_f= input_F_U - yPredicted';    
    copyF_U(input_num+1)=residue_uncertainty_f;
    residuals.uncertainty_f=copyF_U;
    
end