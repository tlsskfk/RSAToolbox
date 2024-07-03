function [ output ] = Compute_Standard_gradCPT_Vars(response)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%Analysis of CPT data


%commissions

commission_errors=size(find(response(:,1)==1 & response(:,7)==-1),1);
correct_omission=size(find(response(:,1)==1 & response(:,7)==0),1);
commission_rate=commission_errors/(correct_omission+commission_errors);


%omissions

omission_errors=size(find(response(:,1)==2 & response(:,7)==0),1);
correct_response=size(find(response(:,1)==2 & response(:,7)==1),1);
omission_rate=omission_errors/(correct_response+omission_errors);

%errors

error_rate=(omission_errors+commission_errors)/(correct_response+omission_errors+commission_errors+correct_omission);


%RT of correct responses

meanRT=mean(response(find(response(:,5)~=0 & response(:,1)==2),5));

%STD of correct responses

STD_RT=std(response(find(response(:,5)~=0 & response(:,1)==2),5));

pHit=1-commission_rate(1,1);
pFA=omission_rate(1,1);
pHit(find(pHit(:)==1))=1-.5/(commission_errors+correct_omission);
pHit(find(pHit(:)==0))=.5/(commission_errors+correct_omission);
pFA(find(pFA(:)==0))=.5/(omission_errors+correct_response);
pFA(find(pFA(:)==1))=1-.5/(omission_errors+correct_response);
HIT=norminv(pHit);
FA=norminv(pFA);
dprime(1,1)=HIT-FA;
criterion(1,1)=(-1*(HIT+FA))/2;
output=[commission_rate,omission_rate,error_rate,dprime,criterion,meanRT,STD_RT];
end

