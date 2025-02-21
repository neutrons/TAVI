clear all;
temp=load('lattice.dat', '-ASCII');
sample.a=temp(1);
sample.b=temp(2);
sample.c=temp(3);
sample.alpha=temp(4)*pi/180;
sample.beta=temp(5)*pi/180;
sample.gamma=temp(6)*pi/180;
lambda=temp(7);

[B,V,Vstar,latticestar]=Bmatrix(sample);

theta20=0.0;
omega0=0.0;
chi0=0.0;

% Read in the observations
data=load('observMorder.dat', '-ASCII');
h=data(:,1);
k=data(:,2);
l=data(:,3);
theta2=(data(:,4)-theta20)*pi/180;
omega=(data(:,5)-data(:,4)/2-omega0)*pi/180;
chi=(data(:,6)-chi0)*pi/180;
phi=data(:,7)*pi/180;

Q=[h(:) k(:) l(:)]';
Qc=B*Q;

% calculate Tc
%h1c=Qc(:,1);
%h2c=Qc(:,2);
%h3c=cross(h1c,h2c);
%t1c=h1c/(sqrt(h1c(1)^2+h1c(2)^2+h1c(3)^2));
%t3c=h3c/(sqrt(h3c(1)^2+h3c(2)^2+h3c(3)^2));
%t2c=cross(t3c,t1c);
%Tc=[t1c(:) t2c(:) t3c(:)];

% Calculate the matrix of u_phi
U(1,:)=cos(omega).*cos(chi).*cos(phi)-sin(omega).*sin(phi);
U(2,:)=cos(omega).*cos(chi).*sin(phi)+sin(omega).*cos(phi);
U(3,:)=cos(omega).*sin(chi);
%h1p=U(:,1);
%h2p=U(:,2);
%ang1=acos(sum(h1p.*h2p))*180/pi;
%h3p=cross(h1p,h2p);
%t1p=h1p/(sqrt(h1p(1)^2+h1p(2)^2+h1p(3)^2));
%t3p=h3p/(sqrt(h3p(1)^2+h3p(2)^2+h3p(3)^2));
%t2p=cross(t3p,t1p);
%Tp=[t1p(:) t2p(:) t3p(:)];

UBmatrix=load('UBmatrix.dat','-ASCII');
%UBmatrix=load('UBmatrix_fmin.dat','-ASCII');

% calculate the (h,k,l) index from observation.
fprintf('Calcuated (H,K,L) from observation:\n');
for i=1:length(theta2)
hphinew=2*sin(theta2(i)/2)/lambda*U(:,i);
H(:,i)=inv(UBmatrix)*hphinew;
fprintf('%5.3f\t%5.3f\t%5.3f\n',H(:,i));
end
