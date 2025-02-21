%writen by F. Ye
warning off
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextFontSize',12)
set(gcf,'PaperOrientation','landscape');


data1=load('scan_MorderAF4_12K.dat','-ASCII');
data2=load('scan_MorderAF4_24K.dat','-ASCII');
foutid=fopen('MnWO4_magneticAF4_12K.int','w');


scannum=data1(:,1);
h=data1(:,2);
k=data1(:,3);
l=data1(:,4);
theta2=data1(:,5);
int1=data1(:,6);
interr1=data1(:,7);
int2=data2(:,6);
interr2=data2(:,7);
% no substracting
%int=int1;
%interr=interr1;
int=int1-int2;
interr=interr1+interr2;
idx=find(int<0);
int(idx)=0;
interr(idx)=0;
width=data1(:,8);
widtherr=data1(:,9);

figure(1); clf;
errorbar(theta2,width,widtherr,'ro');
xlabel('2theta');
ylabel('FWHM');
set(gca,'ylim',[0.2,1.8])

fprintf(foutid,'Single crystal data of MnWO4-13p5Co (hb3a)\n');
%fprintf(foutid,'(3i5,2f8.2,i4,3f8.2)\n');
%fprintf(foutid,'1.53600  0   0\n');
fprintf(foutid,'(4i5,2f8.2,i4,3f8.2)\n');
fprintf(foutid,'1.53600  0   0\n');
fprintf(foutid,'1\n');
fprintf(foutid,'1 0.5 0 0\n');
for i=1:length(h)
    % in shelFORMAT(3I4,2F8.2,I4).  
    %fprintf(foutid,'%4.0f%4.0f%4.0f%8.2f%8.2f%4d\n',h(i),k(i),l(i),int(i),interr(i),0);
    if (width(i)<1.6 &widtherr(i)>0 )
	% generate magnetic int file
	fprintf(foutid,'%5.0f%5.0f%5.0f%5d%8.2f%8.2f%4.0d\n',h(i),k(i),l(i),1,int(i),interr(i),1);
	%fprintf(foutid,'%5.0f%5.0f%5.0f%5d%8.2f%8.2f %5.2f\n',h(i),k(i),l(i),1,int(i),interr(i),theta2(i));
	% generate nuclear int file
	%fprintf(foutid,'%5.0f%5.0f%5.0f%8.2f%8.2f%4.0d\n',h(i),k(i),l(i),int(i),interr(i),1);
     end
 end
fclose(foutid);
fclose all;
