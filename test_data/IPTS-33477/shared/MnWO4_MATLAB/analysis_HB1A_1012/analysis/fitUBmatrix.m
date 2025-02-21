%writen by F. Ye
warning off
clear;

temp=load('UBmatrix.dat', '-ASCII');
UB11=temp(1,1); UB12=temp(1,2); UB13=temp(1,3);
UB21=temp(2,1); UB22=temp(2,2); UB23=temp(2,3);
UB31=temp(3,1); UB32=temp(3,2); UB33=temp(3,3);
theta20=0.00;
chi0= 0.00;

data=load('observ.dat', '-ASCII');
 
%%%%%%%%%%%%%%%Initial Value%%%%%%%%%%%%%%%%%
%  
pa=[UB11,UB12,UB13,UB21,UB22,UB23,UB31,UB32,UB33,theta20,chi0];
ia=[1 1 1 1 1 1 1 1 1 1 1];

[estimates, model] = calcUBfmin(data,pa,ia);
[sse,fval,latticenew,UBmatrix] = model(estimates);

idx=find(ia);
for i=1:length(idx)
    pa(idx(i))=estimates(i);
end
for i=1:length(pa)
    fprintf('pa(%d)=% 5.3e\n',i,pa(i));
end

fval=fval';
for i=1:length(fval(:,1))
fprintf('%5.3f\t%5.3f\t%5.3f\n',fval(i,:));
end

latticenew.alpha=latticenew.alpha*180/pi;
latticenew.beta=latticenew.beta*180/pi;
latticenew.gamma=latticenew.gamma*180/pi;

save UBmatrix_fmin.dat UBmatrix -ASCII
