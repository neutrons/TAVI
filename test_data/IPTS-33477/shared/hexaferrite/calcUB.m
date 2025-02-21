clear all;
temp=load('lattice.dat', '-ASCII');
sample.a=temp(1);
sample.b=temp(2);
sample.c=temp(3);
sample.alpha=temp(4)*pi/180;
sample.beta=temp(5)*pi/180;
sample.gamma=temp(6)*pi/180;
lambda=temp(7);
B=[];
qi=[];

[B,V,Vstar,latticestar]=Bmatrix(sample);

% 4.846	5.721	5.007
% 89.861	91.491	89.736
% Read in the observations
theta20=0.00;
omega0=0.00;
chi0=0.0;
data=load('fit_observ.dat', '-ASCII');
h=data(:,1);
k=data(:,2);
l=data(:,3);
theta2=(data(:,4)-theta20)*pi/180;
omega=(data(:,5)-data(:,4)/2-omega0)*pi/180;
chi=(data(:,6)-chi0)*pi/180;
phi=data(:,7)*pi/180;

Q=[h(:) k(:) l(:)]';
Qc=B*Q;

for i=1:length(h(:))
   x=h(i);
   y=k(i);
   z=l(i);
   qi(i)=sqrt( scalar(x,y,z,x,y,z,latticestar) )/2/pi;
end

% Calculate the matrix of u_phi
U(1,:)=cos(omega).*cos(chi).*cos(phi)-sin(omega).*sin(phi);
U(2,:)=cos(omega).*cos(chi).*sin(phi)+sin(omega).*cos(phi);
U(3,:)=cos(omega).*sin(chi);
Uv1=U(1,:);
Uv2=U(2,:);
Uv3=U(3,:);

% Construct the UB Matrix for calculation
UB(1,1)=sum(h.*h);
UB(1,2)=sum(h.*k);
UB(1,3)=sum(h.*l);
UB(2,1)=UB(1,2);
UB(2,2)=sum(k.*k);
UB(2,3)=sum(k.*l);
UB(3,1)=UB(1,3);
UB(3,2)=UB(2,3);
UB(3,3)=sum(l.*l);

D1(1)=sum(qi.*Uv1.*h(:)');
D1(2)=sum(qi.*Uv1.*k(:)');
D1(3)=sum(qi.*Uv1.*l(:)');
D2(1)=sum(qi.*Uv2.*h(:)');
D2(2)=sum(qi.*Uv2.*k(:)');
D2(3)=sum(qi.*Uv2.*l(:)');
D3(1)=sum(qi.*Uv3.*h(:)');
D3(2)=sum(qi.*Uv3.*k(:)');
D3(3)=sum(qi.*Uv3.*l(:)');

% Solve linear equations for the UB Matrix element
UBmatrix(1,:)=UB\D1(:);
UBmatrix(2,:)=UB\D2(:);
UBmatrix(3,:)=UB\D3(:);

% calculate the (h,k,l) index from observation.
fprintf('Calcuated (H,K,L) from observation:\n');
for i=1:length(theta2)
hphi=2*sin(theta2(i)/2)/lambda*U(:,i);
H(:,i)=inv(UBmatrix)*hphi;
fprintf('%5.3f\t%5.3f\t%5.3f\n',H(:,i));
end

UBtilt=UBmatrix';
InvG=UBtilt*UBmatrix;
G=inv(InvG);
latticenew.a=sqrt(G(1,1));
latticenew.b=sqrt(G(2,2));
latticenew.c=sqrt(G(3,3));
latticenew.alpha=acos(G(2,3)/latticenew.b/latticenew.c)*180/pi;
latticenew.beta=acos(G(1,3)/latticenew.a/latticenew.c)*180/pi;
latticenew.gamma=acos(G(1,2)/latticenew.a/latticenew.b)*180/pi;

%outname=['result.dat'];
%foutid=fopen(outname,'w');
%fprintf(foutid,'#H\tK\tL\tchi_ob\tphi_ob\tchi_cal\tphi_cal\n');
%for i=1:length(h)
        %fprintf(foutid,'%4.2f\t%4.2f\t%4.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\n',...
        %h(i),k(i),l(i),chi(i)*180/pi,phi(i)*180/pi,newchi(i),newphi(i) );
     %end
     %fclose(foutid);

data=load('index.dat', '-ASCII');
%data=load('observ.dat', '-ASCII');
h=data(:,1);
k=data(:,2);
l=data(:,3);
testq=[h(:),k(:),l(:)]';
hphi=UBmatrix*testq;
for i=1:length(h)
   q=sqrt(hphi(1,i).^2+hphi(2,i).^2+hphi(3,i).^2);
   theta(i)=asin(lambda*q/2)*180/pi;
end

fprintf('\nCalculated angle for new requested\n');
fprintf('(  H,  K,  L), 2theta omega chi phi\n');
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
   fprintf('(%5.2f %5.2f %5.2f) mv 2theta %5.2f omega %5.2f chi %5.2f phi %5.2f\n',...
       h(i),k(i),l(i),theta(i)*2+theta20,theta(i)+omega0,newchi(i)+chi0,newphi(i) );
 %
 %  fprintf('(%3d %3d %3d) %5.2f %5.2f %5.2f %5.2f\n',latticenew.a,k(i),l(i),theta(i)*2,theta(i),newchi(i),newphi(i) );
end

save UBmatrix.dat UBmatrix -ASCII
%Add by Huibo to automatically update UBmatrix in UBConf
%UBname=['../UBConf/UBmatrix1.dat'];
%foutubid=fopen(UBname,'w');
%fprintf(foutubid,sprintf('<?xml version="1.0" encoding="utf-8" standalone="no"?>\n'));
%   fprintf(foutubid,sprintf('<!--UB Matrix Information-->\n'));
%   fprintf(foutubid,sprintf('<ubmatrix>\n'));
%   fprintf(foutubid,sprintf('  <time time="10/24/2010 9:43:41 PM" />\n'));
     %lattice_para
% comment=['  <unitcell a="', sprintf('%5.2f',sample.a) '" b="', sprintf('%5.2f',sample.b) '" c="' sprintf('%5.2f',sample.c)...
%         '" alpha="', sprintf('%3.0f',temp(4)) '" beta="', sprintf('%3.0f',temp(5)) '" gamma="', sprintf('%3.0f',temp(6)) '" />'];
 %fprintf(foutubid,'%s\n',comment);
     %wavelenth
%     comment=['  <wavelength lambda="', sprintf('%3.2f',lambda) '" />'];
%     fprintf(foutubid,'%s\n',comment);
   %UBmatrix
comment=['  <matrix matrix="', sprintf('%5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f',UBmatrix(1,:), UBmatrix(2,:), UBmatrix(3,:)) '" />']
%fprintf(foutubid,'%s\n',comment);
%fprintf(foutubid,sprintf('</ubmatrix>\n'));

fprintf('%5.3f\t%5.3f\t%5.3f\n',latticenew.a,latticenew.b,latticenew.c);
fprintf('%5.3f\t%5.3f\t%5.3f\n',latticenew.alpha,latticenew.beta,latticenew.gamma);


%fclose(foutubid);
ubmatrix_line=[sprintf('%5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f',UBmatrix(1,:), UBmatrix(2,:), UBmatrix(3,:))]
