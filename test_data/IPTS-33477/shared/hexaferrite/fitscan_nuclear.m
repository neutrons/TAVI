%writen by F. Ye
warning off
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextFontSize',12)
set(gcf,'PaperOrientation','landscape');
figure(1); clf;

UBmatrix=load('UBmatrix.dat.hexaferrite','-ASCII');
InvUB=inv(UBmatrix);
lambda=2.38;
mcu=9800;
chi0=0;

%finid=fopen('scan_nuclear.index','r');
foutid=fopen('scan_nuclear.dat','w');

%while feof(finid) == 0
for scannum = 999:1211
%for scannum = 1002:1002
%tline = fgets(finid);
  %  scannum=strread(tline,'%d','delimiter',',');
    filename=['/HFIR/HB1A/IPTS-33477/exp1012/Datafiles/HB1A_exp1012_scan' sprintf('%04d',scannum),'.dat'];
    fprintf('%s\n',filename);
	
    [data, xlab, ylab]=spiceloadhb3a(filename);
    scanned_var=std(data(:,1:4));  %find which columns are changing
    %scanned_var=find(scanned_var > 1e-2);  %index those by column number
    scanned_var=find(scanned_var == max(scanned_var));  %index those by column number
    theta2=-data(:,1)/180*pi;
    omega=-(data(:,2)-data(:,1)/2)/180*pi;
    chi=(data(:,3)-chi0)/180*pi;
    phi=data(:,4)/180*pi;
    x=data(:,scanned_var(1));
    monitor=data(:,5);
    tt=sort(monitor);
    monitor=tt(3:end-2);
    monitorave=sum(monitor)/length(monitor);
    ratio=monitorave/mcu;
    y=data(:,6)/ratio;
    err=data(:,7)/ratio;
   
    theta2ave=sum(theta2(:))/length(theta2);
    omegaave=sum(omega(:))/length(omega);
    chiave=sum(chi(:))/length(chi);
    phiave=sum(phi(:))/length(phi);
    U1=cos(omegaave).*cos(chiave).*cos(phiave)-sin(omegaave).*sin(phiave);
    U2=cos(omegaave).*cos(chiave).*sin(phiave)+sin(omegaave).*cos(phiave);
    U3=cos(omegaave).*sin(chiave);
    U=[U1 U2 U3];
    hphiave=2*sin(theta2ave/2)/lambda*U(:);
    H=InvUB*hphiave;
    %have=H(1);
    %kave=H(2);
    %lave=H(3);
    %[have kave lave]
    % AF4 phase
    %have=round((round(H(1)*10)/10-0.236)*2)/2;
    %kave=round((round(H(2)*10)/10-0.5)*2)/2;
    %lave=round((round(H(3)*10)/10+0.438)*2)/2;
     
    % nuclear phase
    have=round(H(1)*10)/10;
    kave=round(H(2)*10)/10;
    lave=round(H(3)*10)/10;
    %[have kave lave];

    %%%%%%%%%%%Draw figures of experimental data %%%%%%%%%%%%%%%%%%%
    figure(1);clf;
    width=0.6;
    height=0.6;
    subplot('Position',[0.20 0.20 width height]);
    errorbar(x,y,err,'ro');
    [bkgrnd,slope,area,center,width]=estwidthhb3(x,y,err);
    pa(1)=bkgrnd; pa(2)=slope; pa(3)=area; pa(4)=center; pa(5)=width;
    slope=0;
    ia=[1 0 1 1 1];
    if area==0
    ia=[0 0 0 0 0];
    end
    [bestpa,bestdpa,chisqN,CN,PQ,nit,kvg,details]=nlfit(x,y,err,'gas_full',pa,ia);
    Eplot=min(x):(max(x)-min(x))/length(x)/4:max(x);
    Fplot=gas_full(Eplot,bestpa);
    %intensity=bestpa(3);
    %interr=bestdpa(3);
    %Lorenz correction
    intensity=bestpa(3)*sin(theta2ave);
    interr=bestdpa(3)*sin(theta2ave);
    width=abs(bestpa(5));
    widtherr=bestdpa(5);
    dx=[x(2:size(x,1))-x(1:(size(x,1)-1))];
    dxave=sum(dx)/length(dx);
    center=bestpa(4)-dxave/2;
    figure(1); hold on; plot(Eplot,Fplot,'k-');
    %titlename=['Scan-#: ' sprintf('%04d (%2.0f,%2.0f,%2.0f)',scannum,have,kave,lave) ];
    titlename=['Scan-#: ' sprintf('%04d (%3.2f,%3.2f,%3.2f)',scannum,have,kave,lave) ];
    title(titlename);
    xlabel('omega (degree)');
    ylabel('counts');
    %pause;
    fprintf(foutid,'%4d%4.0f%4.0f%4.0f  %8.2f %8.2f %8.2f %8.2e %8.2e\n',scannum,have,kave,lave,theta2ave*180/pi,intensity,interr,width,widtherr);
end

fclose(foutid);
%fclose(finid);
