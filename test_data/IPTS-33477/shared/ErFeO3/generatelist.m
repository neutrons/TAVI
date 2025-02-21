
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
srange=1;
maxt=22;
%maxt=6;
mint=2;
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

data=load('scanlist_Xtal.dat', '-ASCII');
h=data(:,1);
k=data(:,2);
l=data(:,3);
int=data(:,4);
int=int(:)';
maxint=max(int);
minint=min(int);
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

idx=find( (theta2>=8 & theta2<89) & ( (newchi > -70)& (newchi < 40))& abs(theta2-61)>3 & abs(theta2-72)>3 );
%output=[qlist(1,idx); qlist(2,idx); qlist(3,idx); theta2(idx); omega(idx); newchi(idx); newphi(idx);]';
output=[qlist(1,idx); qlist(2,idx); qlist(3,idx); theta2(idx); omega(idx); newchi(idx); newphi(idx); int(idx)]';
newoutput=output;
newh=newoutput(:,1);
newk=newoutput(:,2);
newl=newoutput(:,3);
newtheta2=newoutput(:,4);
newomega=newoutput(:,5);
newchi=newoutput(:,6);
newphi=newoutput(:,7);
newint=newoutput(:,8);

data=[newh(:) newk(:) newl(:) newtheta2(:) newomega(:) newchi(:) newphi(:)];
   
% I am generating the list for scanning.
totaltime=0;
outname=['Collecting_HB1A.macro'];
foutid=fopen(outname,'w');
iffirstscan=1;
length(newh)
rate= [45/50,42/40,100/76];
temp = [0,0,0];
for i=1:length(newh)
      comment=['scantitle "Nuclear_@(temp.position)K_2axis: (', sprintf('%4.2f ',newh(i),newk(i),newl(i)), ')"'];
      fprintf(foutid,'%s\n',comment);
      % comment=['br  ', sprintf('%d %d %d',newh(i),newk(i),newl(i))];
      comment = sprintf('mv s2 %.2f s1 %.2f chi %.2f phi %.2f', -newtheta2(i), -newomega(i)+1, newchi(i), newphi(i));
      fprintf(foutid,'%s\n',comment);
      
      diff = [abs(theta2(i)-temp(1)), abs(newchi(i)-temp(2)), abs(newphi(i)-temp(3))];
      time = diff.*rate;
      overhead = max(time)+77;
      roundtime=maxt-round((log10(newint(i))-log10(minint))/(log10(maxint)-log10(minint))*(maxt-mint));
      fprintf(foutid,sprintf('preset time %2d\n',roundtime));
      %pause
      

      %comment = 'scan s1 @(s1)+1.5 @(s1)-1.5 0.15; drive s1 @(com)';
      %fprintf(foutid,'%s\n',comment);
      %comment = 'scan chi @(chi)+5 @(chi)-5 1; drive chi @(com)';
      %fprintf(foutid,'%s\n',comment); 
      %comment = 'th2th 2 -2 0.2';
      comment =  sprintf('scan s1 %.2f %.2f 0.1', -newomega(i)+1, -newomega(i)-1);
      fprintf(foutid,'%s\n',comment);

      
      
      % comment = 'th2th 2 -2 0.2';
      % fprintf(foutid,'%s\n',comment);
      %sroundtime=maxt-round((log10(newint(i))-log10(minint))/(log10(maxint)-log10(minint))*(maxt-mint));
      %fprintf(foutid,sprintf('preset time %2d\n',roundtime));
      %fprintf(foutid,sprintf('scanrel omega %5.2f %5.2f 0.20\n',-srange,srange));
      totaltime=totaltime+roundtime*(round(2*srange/.1)+1)*1.2+overhead;
      
      temp = [theta2(i), newchi(i), newphi(i)];
      
end
fclose(foutid);
fprintf('Total hours needed: %5.2f\n',totaltime/60/60);
