%writen by F. Ye
warning off
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextFontSize',12)
set(gcf,'PaperOrientation','landscape');


data1=load('MnWO4_magneticAF5_5K.int.new','-ASCII');
data2=load('MnWO4_magneticAF5_9K.int.new','-ASCII');
foutid=fopen('MnWO4_magneticAF5_compare.dat','w');

theta5K=data1(:,end);
theta9K=data2(:,end);

[c,ia,ib]=intersect(theta5K,theta9K);
data1new=data1(ia,:);
h1=data1new(:,1);
k1=data1new(:,2);
l1=data1new(:,3);
int1=data1new(:,5);
int1err=data1new(:,6);

data2new=data2(ib,:);
h2=data2new(:,1);
k2=data2new(:,2);
l2=data2new(:,3);
int2=data2new(:,5);
int2err=data2new(:,6);

for i=1:length(c)
    %fprintf(foutid,'%3.0f%3.0f%3.0f%8.2f%8.2f%8.2f%8.2f%3.0f%3.0f%3.0f\n',h1(i),k1(i),l1(i),int1(i),int1err(i),int2(i),int2err(i),h2(i),k2(i),l2(i));
    fprintf(foutid,'%8.2f%8.2f%8.2f%8.2f%3.0f%3.0f%3.0f\n',int1(i),int1err(i),int2(i),int2err(i),h2(i),k2(i),l2(i));
    %fprintf(foutid,'%8.2f%8.2f%8.2f%8.2f\n',int1(i),int1err(i),int2(i),int2err(i));
end
fclose(foutid);
