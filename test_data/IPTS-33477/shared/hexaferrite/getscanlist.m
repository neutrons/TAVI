%general rotation matrix for space group #55, Pbam
% this can be obtained by lookup the space group
a1=[1 0 0; 0 1 0; 0 0 1];
a2=[-1 0 0; 0 -1 0; 0 0 1];
a3=[-1 0 0; 0 1 0; 0 0 -1];
a4=[1 0 0; 0 -1 0; 0 0 -1];
a5=[-1 0 0; 0 -1 0; 0 0 -1];
a6=[1 0 0; 0 1 0; 0 0 -1];
a7=[1 0 0; 0 -1 0; 0 0 1];
a8=[-1 0 0; 0 1 0; 0 0 1];

data=load('GdMn2O5_nuclear.hkl','-ASCII');
% find reflections with intensity greater than 10
%idx=find(data(:,10)<0.1);
idx=find(data(:,10)<200);
%idx=find(data(:,10)<30);
data(idx,:)=[];

qlisttotal=[];
for i=1:length(data(:,1))
    h=data(i,3);
    k=data(i,4);
    l=data(i,5);
    int=data(i,10);
    qlist=[h k l]';
    q1=a1*qlist;
    q2=a2*qlist;
    q3=a3*qlist;
    q4=a4*qlist;
    q5=a5*qlist;
    q6=a6*qlist;
    q7=a7*qlist;
    q8=a8*qlist;
    qlist=[q1 q2 q3 q4 q5 q6 q7 q8]';
    qlisttotal=[qlisttotal; qlist ones(8,1)*data(i,10)];
    end
    qqlist=qlisttotal(:,1:3);
    [q,idx,j]=unique(qqlist,'rows');
    newqlist=qlisttotal(idx,:);

%since lattice parameter is: 7.35, 8.53, 5.68
%first sort agaist k, then agaist h, then l 
%which means varying k first
qlistfinal=sortrows(newqlist,[3 1 2]);

%find only positive h,k, and l
%idx=find( qlistfinal(:,1)>=0 & qlistfinal(:,2)>=0 & qlistfinal(:,3)>=0 );
%idx=find( qlistfinal(:,2)>=0 & qlistfinal(:,3)>=0 );
%idx=find(  qlistfinal(:,2)<=0 & qlistfinal(:,3)>=0 );
%part1
idx=find( qlistfinal(:,3)>=0 );
%idx=find( qlistfinal(:,2)<0 );
data=qlistfinal(idx,:);
h=data(:,1);
k=data(:,2);
l=data(:,3);
int=data(:,4);

foutid=fopen('scanlist_Xtal.dat','w');
for i=1:length(h)
	fprintf(foutid,'%5.0f %5.0f %5.0f %8.1f\n',h(i),k(i),l(i),int(i));
    end
fclose(foutid);
