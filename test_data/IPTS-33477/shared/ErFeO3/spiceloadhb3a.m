function [data,xlab,ylab]=spiceloadhb3a(filename)

% function [data,xlab,ylab]=spiceload(filename)
%
% Load routine for SPICE triple axis data files. 
% This routine either accepts a simple filename input 
% argument, or alternatively a compound file
% specification of the form:
%
% Modified by F.e Ye 07-10-2007

%----- Load data, and attempt to auto-check for new or old tascom files ---------

[data,headertext,headers,defxname,defyname,defxvalue,defyvalue]=spicedata(filename);

%ycol=strmatch(defyname,headers,'exact');
%ycol=strmatch('anode3',headers,'exact');
ycol=strmatch('detector',headers,'exact');
theta2col=strmatch('s2',headers,'exact');
omegacol=strmatch('s1',headers,'exact');
chicol=strmatch('chi',headers,'exact');
phicol=strmatch('phi',headers,'exact');
moncol=strmatch('monitor',headers,'exact');
tempcol=strmatch('tsample',headers,'exact');
%tempcol=strmatch('sample_[b]',headers,'exact');
%timecol=strmatch('time',headers,'exact');

y=data(ycol,:);
theta2=data(theta2col,:);
omega=data(omegacol,:);
chi=data(chicol,:);
phi=data(phicol,:);
monitor=data(moncol,:);
temp=data(tempcol,:);
%time=data(timecol,:);
err=sqrt(y);

y=y(:); err=err(:);
theta2=theta2(:); omega=omega(:);
chi=chi(:); phi=phi(:);
monitor=monitor(:); 
temp=temp(:);
idx=find(err==0);
err(idx)=1;

data=[];
data=[theta2,omega,chi,phi,monitor,y,err,temp];

%----- Create labels
xlab=char(defxname); 
ylab=char(defyname); 
%ymon=[' / ' sprintf('%3.0f secs',mean(time)) '(' num2str(mean(monitor)) ')'];
%ylab=[char(defyname) ymon];

return
