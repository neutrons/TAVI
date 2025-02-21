%writen by F. Ye
warning off
clear all;

mcu=6400;
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextFontSize',12)
set(gcf,'PaperOrientation','landscape');
figure(1); clf;

finid=fopen('Morder.index','r');
foutid=fopen('MorderAF2.dat','w');

while feof(finid) == 0
    tline = fgets(finid);
    scannum=strread(tline,'%d','delimiter',',');
    filename=['../Datafiles/HB3A_exp0107_scan' sprintf('%04d',scannum),'.dat'];
    [data, xlab, ylab]=spiceloadhb3a_T(filename);
    scanned_var=std(data(:,1:4));  %find which columns are changing
    %scanned_var=find(scanned_var > 1e-2);  %index those by column number
    scanned_var=find(scanned_var == max(scanned_var));  %index those by column number
    theta2=data(:,1);
    omega=(data(:,2)-data(:,1)/2);
    chi=data(:,3);
    phi=data(:,4);
    x=omega;
    y=data(:,end-2);
    err=data(:,end-1);
    temp=data(:,end);
    taverage=sum(temp)/length(temp);
    temp=ones(length(x),1)*taverage;
    %if (scannum>517)
       %y=y/6*4;
       %err=err/6*4;
    %end
   
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
    iidex=find(y ~= 0 );
    x=x(iidex);
    y=y(iidex);
    err=err(iidex);
    [bestpa,bestdpa,chisqN,CN,PQ,nit,kvg,details]=nlfit(x,y,err,'gas_full',pa,ia);
    Eplot=min(x):(max(x)-min(x))/length(x)/4:max(x);
    Fplot=gas_full(Eplot,bestpa);
    intensity=bestpa(3);
    interr=bestdpa(3);
    center=bestpa(4);
    centererr=bestdpa(4);
    width=abs(bestpa(5));
    widtherr=bestdpa(5);
    figure(1); hold on; plot(Eplot,Fplot,'k-');
    titlename=['Scan-#: ' sprintf('%04d, T=%3.2f K',scannum(1),taverage) ];
    title(titlename);
    xlabel('omega (degree)');
    ylabel('counts');
    %pause;
    fprintf(foutid,'%4.0f %4.2f %5.2e %5.2e %5.2e %5.2e %5.4e %5.2e\n',scannum(1),taverage,intensity,interr,width,widtherr,center,centererr);
end
fclose(foutid);
fclose(finid);
