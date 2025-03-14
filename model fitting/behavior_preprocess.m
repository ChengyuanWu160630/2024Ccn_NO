function [run,RTs,c,dealdata] = behavior_preprocess(run640)

for i=1:8
    for k=2:5
        for j=1:20
            if run640{i,k}.corr(j) == 0
                run640{i,k}.RT(j) = -1;
                
                if j ~= 20
                    run640{i,k}.RT(j+1) = -1;
                elseif k < 5
                    run640{i,k+1}.RT(1) = -1;
                end
            end
        end
    end
end


congruencytrial=0;
incongruencytrial=0;
for i=1:8
 for k=2:5
     for j=1:20            
     
if run640{i,k}.congruency(j)==0 && run640{i,k}.RT(j)~=-1
    congruencytrial=congruencytrial+1;
    congruencyRT(congruencytrial)=run640{i,k}.RT(j);
end

if run640{i,k}.congruency(j)==1 && run640{i,k}.RT(j)~=-1
    incongruencytrial=incongruencytrial+1;
    incongruencyRT(incongruencytrial)=run640{i,k}.RT(j);
end
     end
 end
end


dealdata.incon_std=std(incongruencyRT);
dealdata.con_std=std(congruencyRT);
dealdata.incon_mean=mean(incongruencyRT);
dealdata.con_mean=mean(congruencyRT);

for i=1:8
 for k=2:5
     for j=1:20  

if run640{i,k}.congruency(j)==0 && run640{i,k}.RT(j)>dealdata.con_mean+2.5*dealdata.con_std && j~=20
    run640{i,k}.RT(j:j+1)=-1; 
end

if run640{i,k}.congruency(j)==1 && run640{i,k}.RT(j)>dealdata.incon_mean+2.5*dealdata.incon_std && j~=20
    run640{i,k}.RT(j:j+1)=-1; 
end       

if run640{i,k}.congruency(j)==0 && k<5 && j==20 && run640{i,k}.RT(j)>dealdata.con_mean+2.5*dealdata.con_std
    run640{i,k}.RT(j)=-1;
    run640{i,k+1}.RT(1)=-1; 
end

if run640{i,k}.congruency(j)==1 && k<5 && j==20 && run640{i,k}.RT(j)>dealdata.incon_mean+2.5*dealdata.incon_std
    run640{i,k}.RT(j)=-1;
    run640{i,k+1}.RT(1)=-1; 
end

if run640{i,k}.congruency(j)==0 && k==5 && j==20 && run640{i,k}.RT(j)>dealdata.con_mean+2.5*dealdata.con_std
    run640{i,k}.RT(j)=-1;
end

if run640{i,k}.congruency(j)==1 && k==5 && j==20 && run640{i,k}.RT(j)>dealdata.incon_mean+2.5*dealdata.incon_std
    run640{i,k}.RT(j)=-1;
end

     end
 end
end
run=run640;

for i=1:8
    trial=0;
 for k=2:5
     for j=1:20 
         trial=trial+1;
         RTs{i}(trial)=run640{i,k}.RT(j);
         c{i}(trial)=run640{i,k}.congruency(j);
     end
 end
end

