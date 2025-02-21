%writen by F. Ye
warning off

%data=load('scan_nuclear_total_sum.dat','-ASCII');
%foutid=fopen('MnWO4_nuclear_5K_total_sum.int','w');
data=load('scan_MorderAF5_5K_part2_sum.dat','-ASCII');
foutid=fopen('magneticAF5_05K/MnWO4_magneticAF5_5K_part2_sum.int','w');

scannum=data(:,1);
h=data(:,2);
k=data(:,3);
l=data(:,4);
int=data(:,5);
interr=data(:,6);

fprintf(foutid,'Single crystal data of MnWO4-13p5Co (hb3a)\n');
%fprintf(foutid,'(3i4,2f8.2,i4,3f8.2)\n');
%fprintf(foutid,'1.536  0   0\n');
fprintf(foutid,'(4i5,2f8.2,i4,3f8.2)\n');
fprintf(foutid,'1.53600  0   0\n');
fprintf(foutid,'1\n');
fprintf(foutid,'1 0.227 0.5 -0.472\n');
for i=1:length(h)
	% generate magnetic int file
	fprintf(foutid,'%5.0f%5.0f%5.0f%5d%8.2f%8.2f%4.0d\n',h(i),k(i),l(i),1,int(i),interr(i),1);
	% generate nuclear int file
	if (int(i)>0)
	%fprintf(foutid,'%5.0f%5.0f%5.0f%8.2f%8.2f%4.0d\n',h(i),k(i),l(i),int(i),interr(i),1);
	%fprintf(foutid,'%4.0f%4.0f%4.0f%8.2f%8.2f%4.0d\n',h(i),k(i),l(i),int(i),interr(i),1);
    end
 end
fclose(foutid);
