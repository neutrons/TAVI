function f=gas_full(x,p,ivar) 
%The form I am using here is identical to the originlab gassian
%form. width obtained is the FWHM, no correction needed.
%
if (nargin <3),ivar=ones(length(p));end %default is vary all 
% p(1) bkgrnd
% p(2) slope
% p(3) area
% p(4) center
% p(5) width
x=x(:); 
bkgrnd=p(1); 
slope=p(2); 
area=p(3);
center=p(4);
width=abs(p(5));
f=bkgrnd+x*slope...
+2*sqrt(log(2)/pi)*area/width*exp(-4*log(2)*((x-center)/width).^2);
