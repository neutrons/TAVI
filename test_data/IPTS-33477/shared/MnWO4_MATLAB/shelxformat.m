%writen by F. Ye
warning off
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextFontSize',12)
set(gcf,'PaperOrientation','landscape');


data=load('scan_nuclear.dat','-ASCII');
%data=load('scan_nuclear_1p003.dat','-ASCII');
%data=load('scan_magnetic_AF5_5K.dat','-ASCII');

foutid=fopen('MnWO4_nuc_300K_2X_HB1A.int','w');
%foutid=fopen('MnWO4_nuclear_5K_1p003.int','w');
%foutid=fopen('MnWO4_magneticAF5_5K.int','w');
%foutid2=fopen('MnWO4_magneticAF5_5K_full.int','w');

%data=sortrows(data,-6);

scannum=data(:,1);
h=data(:,2);
k=data(:,3);
l=data(:,4);
theta2=data(:,5);
int=data(:,6);
interr=data(:,7);
width=data(:,8);
widtherr=data(:,9);

figure(1); clf;
errorbar(theta2,width,widtherr,'ro');
xlabel('2theta');
ylabel('FWHM');
set(gca,'ylim',[0.4,2.0])

fprintf(foutid,'Single crystal data of MnWO4 (hb1a)\n');
fprintf(foutid,'(3i5,2f8.2,i4,3f8.2)\n');
fprintf(foutid,'2,37  0   0\n');
%fprintf(foutid,'(4i5,2f8.2,i4,3f8.2)\n');
%fprintf(foutid,'1.54240  0   0\n');
%fprintf(foutid,'1\n');
%fprintf(foutid,'1 0.236 0.5 -0.438\n');
for i=1:length(scannum)
    % in shelFORMAT(3I4,2F8.2,I4).  
    %fprintf(foutid,'%4.0f%4.0f%4.0f%8.2f%8.2f%4d\n',h(i),k(i),l(i),int(i),interr(i),0);
    if (width(i)<1.2 &widtherr(i)>0 &widtherr(i)<0.2)
	% generate magnetic int file
	%fprintf(foutid,'%5.0f%5.0f%5.0f%5d%8.2f%8.2f%4.0d\n',h(i),k(i),l(i),1,int(i),interr(i),1);
	%fprintf(foutid2,'%5.3f\t%5.3f\t%5.3f\t%5.3f\t%8.2f\t%8.2f\n',h(i)+0.236,k(i)+0.5,l(i)-0.438,theta2(i),int(i),interr(i));
	%fprintf(foutid,'%5.0f%5.0f%5.0f%5d%8.2f%8.2f %5.2f\n',h(i),k(i),l(i),1,int(i),interr(i),theta2(i));
	% generate nuclear int file
	fprintf(foutid,'%5.0f%5.0f%5.0f%8.2f%8.2f%4.0d\n',h(i),k(i),l(i),int(i),interr(i),1);
     end
 end
fclose(foutid);
%fclose(foutid2);
