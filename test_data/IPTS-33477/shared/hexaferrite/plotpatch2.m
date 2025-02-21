%writen by F. Ye
warning off
clear all;

figure; close;
figure(1);clf;
width=0.7;
height=0.7;
subplot('Position',[0.15 0.15 width height]);
set(gcf,'PaperOrientation','portrait');
set(gcf, 'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
width=11.0;	% Initialize a variable for width.
height=8.5;	% Initialize a varible for height.
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

data=load('scanrebin.dat','-ASCII');
x1=data(:,1);
y1=data(:,2);
x2=data(:,3);
y2=data(:,4);
x3=data(:,5);
y3=data(:,6);
x4=data(:,7);
y4=data(:,8);
int=data(:,9);
%idx=find(int==0);
%int(idx)=NaN;

for idx=1:length(x1)-1
    xvec=[x1(idx),x2(idx),x3(idx),x4(idx)];
    yvec=[y1(idx),y2(idx),y3(idx),y4(idx)];
    zvec=[1 1 1 1];
    cvec=[int(idx) int(idx+1) int(idx+1) int(idx)];
    %patch(xvec,yvec,log(cvec+1));
    patch(xvec,yvec,cvec);
    %pause
 end;

xlabel('\omega (^{\circ})');
ylabel('Temperature (K)');
shading flat;
set(gca,'xlim',[-1.5 1.5])
set(gca,'ylim',[4 13.8])
colorbar;
%caxis([50,300]);

fclose('all');
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextFontSize',12);
set(gca,'fontsize',12);
set(get(gca,'xlabel'),'fontsize',12);
set(get(gca,'ylabel'),'fontsize',12);
set(get(gca,'title'),'fontsize',12);

set(gcf, 'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
width=5;   % Initialize a variable for width.
height=7;  % Initialize a varible for height.
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;

myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);
set(gca,'Position',[.12 .15 .7 .7])

%print -deps2 plotCuFeO2_Ga_scanrebin-new.eps
print -dpng plotMnWO4-4p2Co.png
