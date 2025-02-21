%writen by F. Ye
warning off
clear all;
temptotal=[];

set(0,'DefaultAxesFontSize',14);
set(0,'DefaultTextFontSize',14)
set(gcf,'PaperOrientation','portrait');
finid=fopen('Morder.index','r');

figure(1);clf;
width=0.7;
height=0.7;
subplot('Position',[0.15 0.15 width height]);

index=1;
while feof(finid) == 0
    tline = fgets(finid);
    scannum=strread(tline,'%d','delimiter',',');
    filename=['../Datafiles/HB3A_exp0238_scan' sprintf('%04d',scannum),'.dat'];
    fprintf('%s\n',filename);

    [data, xlab, ylab]=spiceloadhb3a_T(filename);
    scanned_var=std(data(:,1:4));  %find which columns are changing
    %scanned_var=find(scanned_var > 1e-2);  %index those by column number
    scanned_var=find(scanned_var == max(scanned_var));  %index those by column number
    theta2=data(:,1);
    omega=(data(:,2)-data(:,1)/2);
    chi=data(:,3);
    phi=data(:,4);
    x=theta2;
    y=data(:,end-2);
    err=data(:,end-1);
    errorbar(x,y,err,'ro');
    temp=data(:,end);
    taverage=sum(temp)/length(temp);
    temp=ones(length(x),1)*taverage;
    temptotal=[temptotal;x(:),temp(:),y(:)];
end
fclose(finid);

[x,y,z]=rebin_grid(temptotal,0.25,0.21);
[c,h]=contourf(x,y,log10(abs(z)+1),50);
%[c,h]=contourf(x,y,z,50);
set(h,'EdgeColor','none');

set(gca,'xlim',[36.9 44.2]);
set(gca,'ylim',[5.3 10.8]);
xlabel('2theta (degree)');
ylabel('Temperature (K)');

set(gcf, 'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
height=7;  % Initialize a varible for height.
width=5;   % Initialize a variable for width.
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;

myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);
set(gca,'Position',[.15 .15 .7 .7])
%grid
%grid minor

print -dpng Morder_10pCu.png
