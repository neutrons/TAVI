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

UU=0.2942;
VV=-0.1724;
WW=0.2146;
GG=-0.1365;

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
   for h=0.1:0.01:0.4;
%   for l=-5.2:0.02:3.2;
       qlist=[qlist;h,1.5,-2*h];
       %qlist=[qlist;1,0,1+l];
   end
h=qlist(:,1);
k=qlist(:,2);
l=qlist(:,3);

%data=load('highInt.index', '-ASCII');
%h=data(:,1);
%k=data(:,2);
%l=data(:,3);

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

%idx=find(theta2>5 & theta2<97.8 & newchi>0 );
%idx=find(theta2>5 & theta2<60 );
idx=find(theta2>5 & theta2<88 & omega>2);
%idx=find( (theta2>5 & theta2<80) & (newphi>-60 & newphi<120) &(newchi<68 & newchi >-20) );
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
outname=['qscanh03h.macro'];
foutid=fopen(outname,'w');
iffirstscan=1;
comment=['scantitle "MnWO4 2p-Co T@(tsample.position)1.01_mask20_50s_mb=@(mbend)"'];
fprintf(foutid,'%s\n',comment);
comment=['scanon'];
fprintf(foutid,'%s\n',comment);
for i=1:length(newh)
     % tanx=tan(newomega(i)/180*pi);
      %cosx=cos(newomega(i)/180*pi);
      %f=UU.*tanx.*tanx+VV.*tanx+WW+GG./cosx./cosx;
      %range=sqrt(f)*12;
      % for perfect sample.
      % range=sqrt(f)*4;
      
     % comment=[sprintf('%4.2f ',newh(i),newk(i),newl(i))];
      %fprintf(foutid,'%s',comment);
      comment=['drive 2theta ', sprintf('%5.2f',newtheta2(i)) ...
      ' omega ', sprintf('%5.2f',newomega(i))...
      ' chi ', sprintf('%5.2f',newchi(i)) ' phi ',sprintf('%5.2f',newphi(i))];
      fprintf(foutid,'%s\n',comment);
      comment=['count 2'];
      fprintf(foutid,'%s\n',comment);
      comment=['stepend'];
      fprintf(foutid,'%s\n',comment);
      
      oldpos=[newtheta2(i),newomega(i),newchi(i),newphi(i)];
      
end
comment=['scanoff'];
fprintf(foutid,'%s\n',comment);
fclose(foutid);
outname=['qlistqscanh03h.dat'];
fid=fopen(outname,'w');
for i=1:length(qlist(1,:))
comment=[sprintf('%5.2f ',qlist(1,i)),sprintf('%5.2f ',qlist(2,i)),sprintf('%5.2f ',qlist(3,i))];
fprintf(fid,'%s\n',comment);
end
fclose(fid);
