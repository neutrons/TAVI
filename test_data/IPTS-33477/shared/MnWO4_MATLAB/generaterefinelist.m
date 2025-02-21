clear all;
% load lattice information
temp=load('lattice.dat', '-ASCII');
sample.a=temp(1);
sample.b=temp(2);
sample.c=temp(3);
sample.alpha=temp(4)*pi/180;
sample.beta=temp(5)*pi/180;
sample.gamma=temp(6)*pi/180;
lambda=temp(7);
[B,V,Vstar,latticestar]=Bmatrix(sample);
%offset
theta20=0.00;
omega0=0.0;
chi0=0.0;
srange=2.40;
maxt=10;
mint=10;
UBmatrix=load('UBmatrix.dat','-ASCII');

qx=[1,0,0];
qy=[0,1,0];
qz=[0,0,1];
Newq=[qx(:),qy(:),qz(:)]';
Qc=B*Newq;
Qindex=[];

for i=1:3
   q(i)=sqrt(Qc(1,i)^2+Qc(2,i)^2+Qc(3,i)^2);
   thetamax=45*pi/180;
   Nmax=2*sin(thetamax)/lambda/q(i);
   Qindex=[Qindex; ceil(Nmax)];
end
Qindex
   
data=load('scan_nuclear.dat', '-ASCII');
h=data(:,2);
k=data(:,3);
l=data(:,4);

threshold = 60;
key_set = zeros(0, 3);

% Add keys until the threshold is reach
idx = 1;
while size(key_set, 1) < threshold
    new_key = [h(idx), k(idx), l(idx)]; 
    if ~ismember(new_key, key_set, 'rows')
        key_set = [key_set; new_key]; 
    end
    new_key = [-h(idx), -k(idx), -l(idx)]; 
    if ~ismember(new_key, key_set, 'rows')
        key_set = [key_set; new_key]; 
    end
    new_key = [h(idx), -k(idx), l(idx)]; 
    if ~ismember(new_key, key_set, 'rows')
        key_set = [key_set; new_key]; 
    end
    new_key = [-h(idx), k(idx), -l(idx)]; 
    if ~ismember(new_key, key_set, 'rows')
        key_set = [key_set; new_key]; 
    end
    idx = idx+1;
end

h = key_set(:, 1);
k = key_set(:, 2);
l = key_set(:, 3);

% maxint=max(int);
% minint=min(int);
qlist=[h(:) k(:) l(:)];

qlist=qlist';
Qc=B*qlist;
for i=1:length(h)
q(i)=sqrt(Qc(1,i)^2+Qc(2,i)^2+Qc(3,i)^2);
end
theta=asin(lambda*q/2);
theta2=theta*2*180/pi;
omega=theta2/2;

hphi=UBmatrix*qlist;
for i=1:length(h(:))
   newphi(i)=atan(hphi(2,i)/hphi(1,i))*180/pi;
   if newphi(i)>0
      if hphi(1,i)<0
         newphi(i)=newphi(i)-180;
      end
   else
      if hphi(1,i)<0
         newphi(i)=newphi(i)+180;
      end
   end
   newchi(i)=atan(hphi(3,i)/sqrt(hphi(1,i)^2+hphi(2,i)^2))*180/pi;
end

%idx=find( (theta2>=-2 & theta2<90) & ( (newchi > -70)& (newchi < 40)));
idx=find( (theta2>=-2 & theta2<90) & ( (newchi > -70)& (newchi < 40))& (abs(theta2-62)>3 & abs(theta2+62)>3 & abs(theta2-71)>3 & abs(theta2+71)>3 ) );
output=[qlist(1,idx); qlist(2,idx); qlist(3,idx); theta2(idx); omega(idx); newchi(idx); newphi(idx);]';
newoutput=output;
newh=newoutput(:,1);
newk=newoutput(:,2);
newl=newoutput(:,3);
newtheta2=newoutput(:,4);
newomega=newoutput(:,5);
newchi=newoutput(:,6);
newphi=newoutput(:,7);

data=[newh(:) newk(:) newl(:) newtheta2(:) newomega(:) newchi(:) newphi(:)];
   
% I am generating the list for scanning.
totaltime=0;
outname=['Optimizing_HB1A.macro'];
foutid=fopen(outname,'w');
iffirstscan=1;
length(newh)
for i=1:length(newh)
      comment=['scantitle "Nuclear_@(temp.position) (', sprintf('%4.2f ',newh(i),newk(i),newl(i)), ')"'];
      fprintf(foutid,'%s\n',comment);
      % comment=['br  ', sprintf('%d %d %d',newh(i),newk(i),newl(i))];
      comment = sprintf('mv s2 %.2f s1 %.2f chi %.2f phi %.2f', -newtheta2(i), -newomega(i), newchi(i), newphi(i));
      fprintf(foutid,'%s\n',comment);
      comment = 'scan s1 @(s1)+1.5 @(s1)-1.5 0.15; drive s1 @(com)';
      fprintf(foutid,'%s\n',comment);
      comment = 'scan chi @(chi)+5 @(chi)-5 1; drive chi @(com)';
      fprintf(foutid,'%s\n',comment); 
      comment = 'th2th 2 -2 0.2';
      fprintf(foutid,'%s\n',comment);
      % comment = 'th2th 2 -2 0.2';
      % fprintf(foutid,'%s\n',comment);
      % sroundtime=maxt-round((log10(newint(i))-log10(minint))/(log10(maxint)-log10(minint))*(maxt-mint));
      % fprintf(foutid,sprintf('preset time %2d\n',roundtime));
      % fprintf(foutid,sprintf('scanrel omega %5.2f %5.2f 0.20\n',-srange,srange));
      % totaltime=totaltime+roundtime*(round(2*srange/.20)+1)*1.2;
      
end
fclose(foutid);
