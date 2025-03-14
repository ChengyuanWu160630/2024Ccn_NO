function[residuals]=calculate_residue(pwe,alpha,f,congruency,uncertainty,need_to_changeZ)
pwe_noneed_deleted=find(need_to_changeZ==1);
input_alpha=[];
input_f=[];
input_pwe=[];
input_num=[];
input_congruency=[];
input_uncertainty=[];
for i=1:length(pwe_noneed_deleted)
    extract_num=pwe_noneed_deleted(i);
    if isempty(find([80:80:640]==extract_num, 1)) && ~isempty(find(pwe_noneed_deleted==(extract_num+1), 1))
        input_pwe(end+1)=pwe(extract_num);
        input_congruency(end+1)=congruency(extract_num);
        input_f(end+1)=f(extract_num+1);
        input_alpha(end+1)=alpha(extract_num+1);
        input_uncertainty(end+1)=uncertainty(extract_num+1);
        input_num(end+1)=extract_num;
    end
end
   %% PE
    X = [ones(size(input_alpha))', input_alpha', input_f', input_congruency',(input_alpha.*input_f)',(input_alpha.*input_congruency)',(input_alpha.*input_pwe)'];
    b = regress(input_pwe', X);
    yPredicted = X * b;
    residue_pwe= input_pwe - yPredicted';
    pwe(input_num)=residue_pwe;
    residuals.pwe =pwe;   
   %% f
    X = [ones(size(input_alpha))', input_alpha', input_pwe', input_congruency',(input_alpha.*input_f)',(input_alpha.*input_congruency)',(input_alpha.*input_pwe)'];
    b = regress(input_f', X);
    yPredicted = X * b;
    residue_f= input_f - yPredicted';
    f(input_num+1)=residue_f;
    residuals.f =f;   
   %% alpha
    X = [ones(size(input_alpha))', input_pwe', input_f', input_congruency',(input_alpha.*input_f)',(input_alpha.*input_congruency)',(input_alpha.*input_pwe)'];
    b = regress(input_alpha', X);
    yPredicted = X * b;
    residue_alpha= input_alpha - yPredicted';    
    alpha(input_num+1)=residue_alpha;
    residuals.alpha=alpha;
    residuals.input_num=input_num;
   %% alpha_pwe
    X = [ones(size(input_alpha))', input_alpha',input_pwe', input_f', input_congruency',(input_alpha.*input_f)',(input_alpha.*input_congruency)'];
    alpha_pwe=1:length(alpha);
    b = regress((input_alpha.*input_pwe)', X);
    yPredicted = X * b;
    residue_alpha_pwe= (input_alpha.*input_pwe)- yPredicted';    
    alpha_pwe(input_num+1)=residue_alpha_pwe;
    residuals.alpha_pwe=alpha_pwe;   
   %% uncertainty
    X = [ones(size(input_alpha))',input_alpha', input_pwe', input_f', input_congruency',(input_alpha.*input_f)',(input_alpha.*input_congruency)',(input_alpha.*input_pwe)'];
    b = regress(input_uncertainty', X);
    yPredicted = X * b;
    residue_uncertainty= input_uncertainty - yPredicted';    
    uncertainty(input_num+1)=residue_uncertainty;
    residuals.uncertainty=uncertainty;   
end