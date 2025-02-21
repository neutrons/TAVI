%writen by F. Ye
warning off

UBmatrix=load('UBmatrix.dat','-ASCII');
InvUB=inv(UBmatrix);
lambda=1.536;

theta20=0.42;
omega0=0.0;
chi0=-0.11;

scanlist=load('nuclear_scanlist.dat','-ASCII');
%scanlist=load('magnetic_scanlist.dat','-ASCII');
ii=1;

finid=fopen('scan_nuclear.index','r');
%finid=fopen('scan_Morder5K.index','r');

foutid=fopen('scan_nuclear_list.dat','w');
%foutid=fopen('scan_Morder5K_list.dat','w');

%foutid2=fopen('scan_nuclear_goodlist.dat','w');
foutid2=fopen('scan_nuclear_badlist.dat','w');
%foutid2=fopen('scan_magnetic_badlist.dat','w');

while feof(finid) == 0
    tline = fgets(finid);
    scannum=strread(tline,'%d','delimiter',',');
    filename=['../Datafiles/HB3A_exp0107_scan' sprintf('%04d',scannum),'.dat'];
    %fprintf('%s\n',filename);
	
    [data, xlab, ylab]=spiceloadhb3a(filename);
    scanned_var=std(data(:,1:4));  %find which columns are changing
    %scanned_var=find(scanned_var > 1e-2);  %index those by column number
    scanned_var=find(scanned_var == max(scanned_var));  %index those by column number
    theta2=data(:,1)/180*pi;
    omega=(data(:,2)-data(:,1)/2)/180*pi;
    chi=data(:,3)/180*pi;
    phi=data(:,4)/180*pi;
    x=data(:,scanned_var(1));
    y=data(:,end-1);
    err=data(:,end);
   
    theta2ave=sum(theta2(:))/length(theta2)+theta20*pi/180;
    omegaave=sum(omega(:))/length(omega);
    chiave=sum(chi(:))/length(chi)+chi0*pi/180;
    phiave=sum(phi(:))/length(phi);
    U1=cos(omegaave).*cos(chiave).*cos(phiave)-sin(omegaave).*sin(phiave);
    U2=cos(omegaave).*cos(chiave).*sin(phiave)+sin(omegaave).*cos(phiave);
    U3=cos(omegaave).*sin(chiave);
    U=[U1 U2 U3];
    hphiave=2*sin(theta2ave/2)/lambda*U(:);
    H=InvUB*hphiave;

    % nuclear phase
    have=H(1);
    kave=H(2);
    lave=H(3);
    % magnetic phase
    %have=H(1)-0.214;
    %kave=H(2)-0.5;
    %lave=H(3)+0.457;

    HH=scanlist(ii,1);
    KK=scanlist(ii,2);
    LL=scanlist(ii,3);
    err=(H(1)-HH)^2+(H(2)-KK)^2+(H(3)-LL)^2;
    fprintf(foutid,'%4d\t%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t%6.3f\n',scannum,H(1),H(2),H(3),HH,KK,LL,err);
    if (err>0.01)
    fprintf(foutid2,'%4d\t%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t%6.3f\n',scannum,H(1),H(2),H(3),HH,KK,LL,err);
    end
    ii=ii+1;
end

fclose(foutid);
fclose(foutid2);
fclose(finid);
