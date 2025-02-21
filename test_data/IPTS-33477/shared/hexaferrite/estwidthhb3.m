function [bkgrnd,slope,area,center,width]=estwidthhb3(x,y,err)
% [bkgrnd,slope,area,center,width]=estwidth(x,y,err)
%estimate peak center, peak width and peak height

x=x(:);
y=y(:);
err=err(:);
data=[x,y,err];
newdata=sortrows(data,1);
x=newdata(:,1);
y=newdata(:,2);
err=newdata(:,3);

data1=[x(:) y(:)];
newdata=sortrows(data1,2);
y1=newdata(:,2);
y1=sortrows(y1(:));
idx1=round(length(y)/11);
idx2=round(length(y)/4);
newx=newdata(idx1:idx2,1);
newy=newdata(idx1:idx2,2);
[pin1,ny] = polyfit(newx,newy,0);

%bkgrnd=pin1(1);
bkgrnd=sum(newy)/length(newy);
%pause
slope=0;
y=y-(slope*x+bkgrnd);

dx=[x(2)-x(1);x(3:size(x,1))-x(1:(size(x,1)-2));x(size(x,1))-x(size(x,1)-1)];
totarea=dx'*y;
maxx=max(x); minx=min(x);
av = totarea/(maxx - minx);
xcom = sum(x)/length(x);
%for standard gaussian, av=0.295*peak height
x2=x(3:end-3);
y2=y(3:end-3);
center=sum(x2.*y2)/sum(y2);
moment2=sum((x2-center).*(x2-center).*y2)/sum(y2);
sigma=sqrt(moment2)*1.50;

idx=find(x<=center);
% if peak is the first point
ypeak=y(length(idx));
if (length(idx)>1 & length(idx)<length(x))
ypeak=(y(length(idx)-1)+y(length(idx))+y(length(idx)+1))/3;
end
peak=ypeak-av;

%for standard gaussian, peak=0.705*peak height
newidx=find(y<av);
peakarea=totarea-dx(newidx)'*y(newidx)*0.5;
area=peakarea*1.1; %integral from -FWHM/2 to FWHM/2 is 90.26% of total weight
%width = 0.66*area/abs(peak);

width=sigma
y3=sortrows(y(:));
newypeak=sum(y3(end-2:end))/3;

%center
%width
%area

% I am printing out the debug info.
stdbk=std(y3(round(length(y)*0.10):round(length(y)*0.35)));
range=(max(x)-min(x))/4.5;
centeridx=find(abs(x-center)<0.20*range);
peakheight=sum(y(centeridx))/length(centeridx);
idx1=find(abs(x-center) <= range);
idx2=find(abs(x-center) >= range);
sigy=y(idx1);
bky=y(idx2);
std(sigy)/std(bky)
peakheight/stdbk

%if(abs(std(sigy)/std(bky))>0.5 & peakheight/stdbk>4 & (abs(center-xcom)<0.3) )
if( peakheight/stdbk>4 & (abs(center-xcom)<0.3) )
else
  area=0;
  width=2.5;
end

