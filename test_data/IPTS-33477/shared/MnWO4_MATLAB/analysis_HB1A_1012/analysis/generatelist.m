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
maxt=9;
mint=3;
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
   
% This is magnetic Bragg generation.
qlist=[];
   for  h=-9:9
      for k=-6:-1
         for l=0:10
             %if h*l~=0
             %if h==0 & l~=0 & mod(abs(k+l),2)==0
             %if l==0 & mod(abs(k+h),2)==0
             %if mod(abs(-h+k+l),3)==0
            %Spc C_2/c No15
            %if mod(h+k,2)==0 & k~=0 || k==0 & mod(h,2)==0 & mod(l,2)==0;
            %I4/mmm 
            %if mod(h+k+l,2)==0;
                qlist=[qlist;h,k,l];
            %end
         end
      end
   end
h=qlist(:,1);
k=qlist(:,2);
l=qlist(:,3);

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

idx=find( (theta2>=10 & theta2<90) & ( (newchi > -65)& (newchi < 65))& ( (omega-srange)>-41 & (omega+srange)<41 ) );
output=[qlist(1,idx); qlist(2,idx); qlist(3,idx); theta2(idx); omega(idx); newchi(idx); newphi(idx); int(idx)]';
%output=[qlist(1,idx); qlist(2,idx); qlist(3,idx); theta2(idx); omega(idx); newchi(idx); newphi(idx);]';
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
outname=['Collecting_4K_1p542_part2.macro'];
foutid=fopen(outname,'w');
iffirstscan=1;
for i=1:length(newh)
      comment=['scantitle "Nuclear_@(tsample.position) (', sprintf('%4.2f ',newh(i),newk(i),newl(i)), ')"'];
      fprintf(foutid,'%s\n',comment);
      comment=['br  ', sprintf('%d %d %d',newh(i),newk(i),newl(i))];
      fprintf(foutid,'%s\n',comment);
      roundtime=maxt-round((log10(newint(i))-log10(minint))/(log10(maxint)-log10(minint))*(maxt-mint));
      fprintf(foutid,sprintf('preset time %2d\n',roundtime));
      fprintf(foutid,sprintf('scanrel omega %5.2f %5.2f 0.20\n',-srange,srange));
      totaltime=totaltime+roundtime*(round(2*srange/.20)+1)*1.2;
      
end
fclose(foutid);
fprintf('Total hours needed: %5.2f\n',totaltime/60/60);
