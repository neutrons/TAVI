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

% load UBmatrix
UBmatrix=load('UBmatrix.dat','-ASCII');

theta20=0.0;
omega0=0.0;
chi0=-0.0;
srange=2.40;

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
   
% This is magnetic Bragg generation.
qlist=[];
   for  l=0:6
      for h=0:7
         for k=-8:8
             qlist=[qlist;h+0.5,k,l];
         end
      end
   end
h=qlist(:,1);
k=qlist(:,2);
l=qlist(:,3);

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

idx=find( theta2>=9 & theta2<60 &(newchi<50 & newchi >-50) & ( (omega-srange)>-41 & (omega+srange)<41 ) );
output=[qlist(1,idx); qlist(2,idx); qlist(3,idx); theta2(idx); omega(idx); newchi(idx); newphi(idx)]';
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
outname=['Collecting_magnetic_CM_5_55.macro'];

foutid=fopen(outname,'w');
fprintf(foutid,sprintf('preset time 12\n'));
for i=1:length(newh)
   comment=['scantitle "Magnetic @(tsample.position) (', sprintf('%5.3f ',newh(i),newk(i),newl(i)), ')"'];
   fprintf(foutid,'%s\n',comment);
   comment=['br  ', sprintf('%5.1f %5.1f %5.1f',newh(i),newk(i),newl(i))];
   fprintf(foutid,'%s\n',comment);
   fprintf(foutid,sprintf('scanrel omega %5.2f %5.2f 0.20\n',-srange,srange));
   totaltime=totaltime+12*(round(2*srange/.20)+1)*1.2;
      
end
fclose(foutid);
fprintf('Total hours needed: %5.2f\n',totaltime/60/60);
