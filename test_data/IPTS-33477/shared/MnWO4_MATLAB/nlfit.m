function [pa,dpa,chisqN,CN,PQ,nit,kvg,details]=nlfit(x,y,s,func,pa,ia,nitmax,tol,udiff,dtol) 
%[pa,dpa,chi2N,CN,PQ,nit,kvg,details]=nlfit(x,y,s,'func',pa,ia,nitmax,tol,udiff,dtol) 
% 
% Levenberg-Marquardt least squares fit to nonlinear function 
% see NR 15.5, also function Lfit.m 
% Subfunctions used: marqit (calculates the chisquared and needed matrices) 
%                    dfdp (for numerical calculation of derivatives) 
%                    'func' (user supplied fitting function) - see notes below 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%input: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%compulsory input parameters: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% x - independent variable, (may be vector or multidimensional column vector form depending on func) 
% y - measured value, and s, error (vectors each of equal length) 
% func - function name (see details below) 
% pa - guess of initial values for parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%optional input parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ia - list with zero elements corresponding to fixed parameters  DEFAULT: all varied ia=ones 
% nitmax - the maximum number of iterations allowed before the program exits (default 20) 
% tol - the convergence criterion: convergence declared when chi-squared decreases by a 
%       relative amount less than tol (default tol=0.001) 
% udiff - udiff=1 user func supplies derivatives.  udiff=0 (or ~1 ) numerical derivatives 
%         default is udiff =0 (numerical derivatives) 
% dtol - parameter for relative change in numerical derivatives, default dtol=1e-5 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%output: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% pa,dpa - (parameters and errors) 
% chi2N - (chi-squared normalized by degrees of freedom) 
% CN  - normalized correlation matrix - includes varying parameters only (in order) 
% note : CN(i,j)=C(i,j)/sqrt(C(i,i)*C(j,j)) 
%        dpa(i)=sqrt(c(i,i)) 
% PQ - probability that chi2N exceeds that observed (see NR chapter 15). 
% nit - number of iterations to convergence 
% kvg - convergence flag: 
%       0  did not converge due to too many iterations nit >= nitmax 
%       1  converged normally 
%       2  questionable convergence (final lamda > 1e-3 - may indicate bad starting point for pa) 
%  
% details - optional output structure with fields: 
%       chisq (raw chi-squared) 
%       DF (degrees of freedom) 
%       Ndata (number of data points) 
%       Npar (total number of parameters) 
%       Nvar  (number of parameters varied) 
%       DF (degrees of freedom) 
%       C (raw covariance matrix) 
%       final_lamda (final value of lamda used in the fitting process) 
%       yf (fitted values of function at input points) 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%function details: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%user supplied derivatives: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%[f,dfp]=func(x,p,ivar) 
% input: independent variable x, parameters p 
% optional input - list ivar - if omitted the default is ivar = ones(length(p)) 
% in all cases the function outputs f (should be as a column vector) 
% if dfp is requested func calculates the derivatives as well 
% these should be returned in column matrix form with row labels x and column labels p 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%numerical derivatives: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% use f=func(x,p,ivar) to output f 
% The subfunction dfdp calculates the derivatives numerically 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% under development by S.E. Nagler, ORNL, 1999 
% 
 
% check defaults and initial conditions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
y=y(:);%ensure that y is a column vector 
if (isempty(s)),s=ones(size(y));else,s=s(:);end %also ensure that s is a column vector 
% ia is a list with zero elements corresponding to the fixed parameters. 
if (nargin < 6 | isempty(ia)),ia=ones(1,length(pa)); end %default is to vary everything 
if (nargin < 7 | isempty(nitmax)), nitmax=1000; end %default limit on iterations 
if (nargin < 8 | isempty(tol)),tol=0.001;,end %default tolerance determining convergence 
if (nargin < 9 | isempty(udiff)),udiff=0;end %default is numerical derivatives 
if (nargin < 10 | isempty(dtol)), dtol=1e-5;,end %default derivative tolerance 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
dpa=zeros(size(ia)); %default error is zero (applies to fixed parameters) 
ivar=find(ia); %list of elements corresponding to varying parameters 
DF=length(y)-length(ivar); % degrees of freedom 
lamda = 0.003; % initial Marquardt parameter 
nit=0; % initial number of iterations 
chisq_old=realmax; %initial chisq_old huge to ensure iteration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% function MARQIT: calculate alpha, beta, chi_squared and new marquardt parameter alpha 
% beta(k)=-sum_i {(y-f)*dfp(k)/s^2} 
% alpha(k,l) = sum_i {dfp(k).*dfp(l)./s.^2} 
[chisq,alpha,beta]=marqit(func,x,y,s,pa,ivar,udiff,dtol);% initial call 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% begin iterations 
% debugging
while (abs(chisq_old-chisq)> 0.001*chisq) 
   if (nit >=nitmax),break,end %punt if too many iterations 
% solve the equations [alpha*dp=beta modified by multiplying diag(alpha) by 1+lamda 
% dpa(ivar)=inv(alpha')*beta but alpha'\beta is more accurate 
% alpha'(i,i)=alpha(i,i)*(1+lamda) 
		dpa(ivar)=(alpha+lamda*diag(diag(alpha)))\beta(:);  
		pt=pa+dpa; % new trial parameters 
   	[chisq_trial]=marqit(func,x,y,s,pt,ivar,udiff,dtol);% get trial chi_squared 
   if (chisq_trial > chisq) % chisq increases ? 
      lamda=10*lamda; % increase lamda and try again 
      if(lamda > 1e13),chisq_old=chisq;end %punt if lamda is stuck beyond reason 
   else %chisq decreased 
     	chisq_old=chisq; %update old chi-squared 
     	lamda=lamda/10; %decrease lamda 
     	pa=pt; % update parameters 
     	nit=nit+1; % call it an iteration 
     	[chisq,alpha,beta]=marqit(func,x,y,s,pa,ivar,udiff,dtol); %recalculate alpha,beta,chi-squared 
     	disp([nit chisq/DF lamda]) 
  end     
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% check convergence: 
if (nit >= nitmax | lamda > 1e13),kvg=0; 
elseif (lamda > 0.001),kvg=2; 
else kvg=1; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% now that converged or kicked out calculate error quantities with lamda=0 
chisqN=chisq/DF; %normalized chi-squared
C=inv(alpha); % raw correlation matrix
dpa(ivar)=sqrt(chisqN*diag(C)); % error in a(j) is sqrt(Cjj)*sqrt(chisqN)
CN=C./sqrt(abs(diag(C)*diag(C)')); % normalized correlation matrix C(ij)/sqrt[C(ii)*C(jj)]
% PQ is probability that chi_2 exceeds that observed - see NR 6.2.3
% note matlab's definition of gammainc switches arguments compared to NR
% calculate this only if requested since it takes a little time
if (nargout >=5),PQ=1-gammainc(chisq/2,DF/2), end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% more details of the fit if requested 
if (nargout==8), 
	details = struct('chisq',1,'Ndata',1,'Npar',1,'Nvar',1,'DF',1,'C',1,'final_lamda',1,'yf',zeros(size(y))); 
	details.chisq=chisq; %raw chi-squared 
	details.Ndata=length(y); % number of data points 
	details.Npar=length(pa); %total number of parameters 
	details.Nvar=length(ivar); %number of parameters varied 
	details.DF=DF; %degrees of freedom 
	details.C=C; %raw covariance matrix 
	details.final_lamda=lamda; %final value of lamda used in the fitting process 
	details.yf=feval(func,x,pa,ivar); %fitted values of the function at input points 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
return 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
function [chisq,alpha,beta] = marqit(func,x,y,s,pa,ivar,udiff,dtol) 
% beta(k)=-sum_i {(y-f)*dfp(k)/s^2} 
% alpha(k,l) = sum_i {dfp(k).*dfp(l)./s.^2} 
switch nargout % what is requested: chisq or full treatment? 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
case 1 %only chisq requested 
   f=feval(func,x,pa,ivar); 
   wdiff=(y(:)-f(:))./s(:); %weighted difference 
   chisq=sum(wdiff.^2); % chi_squared 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
otherwise %get alpha and beta as well 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   switch udiff % user supplied derivatives? 
   case 1 % derivatives supplied by user 
      [f,dfp]=feval(func,x,pa,ivar); % calculate functions and derivatives 
   otherwise %calculate derivatives numerically 
      f=feval(func,x,pa,ivar); 
      dfp=dfdp(func,x,pa,ivar,dtol); 
   end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   f=f(:); %ensure that f is a column vector 
   wdiff=(y(:)-f(:))./s(:); %weighted difference 
   chisq=sum(wdiff.^2); % chi_squared 
   NP=length(ivar);%number of varying parameters 
   % beta(k)=sum_i {(y-f)*dfp(k)/s^2} (column vector) 
   beta=sum(repmat(wdiff./s(:),1,NP).*dfp)'; 
   % set up alpha(k,l) = sum_i {dfp(k).*dfp(l)./s.^2} 
   %diagonal is multiplied by (1+lamda) 
   alpha=dfp./repmat(s(:),1,NP); %normalize derivative by sigma 
   alpha=alpha'*alpha; %product alpha'*alpha to get the correct sum 
   % alpha=alpha+lamda*diag(diag(alpha));%add lamda times diagonal in calling function 
end 
return 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
function [dfp]=dfdp(func,x,pa,ivar,dtol) 
%calculate derivatives using a two-sided numerical method 
if(nargin < 5),dtol = 1e-8;end % default tolerance for derivatives 
[np,nx]=size(x); %just in case x multi-valued 
if(np==1 | nx==1),x=x(:);np=length(x);end %ensure that x is a column vector if need be 
dfp=zeros(np,length(ivar)); %preallocate space for dfdp as zeros 
for n=1:length(ivar) %loop through varying parameters 
   h=zeros(size(pa)); %initialize h 
   % make difference h an exactly representable number 
   t=pa(ivar(n))+dtol*pa(ivar(n)); 
   h(ivar(n))=t-pa(ivar(n));  
   if(pa(ivar(n))==0),h(ivar(n))=1+1e-8-1;end %protect against zero values 
   % calculate the derivative 
   dfp(:,n)=(feval(func,x,pa+h,ivar)-feval(func,x,pa-h,ivar))/(2*h(ivar(n))); 
end 
return 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   end of core functions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [f,dfp]=gasdum(x,p,ivar) 
%f=bg+A*exp[-(x-a)^2/(2*s^2)]/sqrt(2*pi*s) 
if (nargin <3),ivar=ones(length(p));end %default is vary all 
x=x(:); 
A=p(1); 
a=p(2); 
s=abs(p(3)); 
b0=p(4); 
z=(x-a)./s; 
f=b0+A*exp(-z.^2/2)/sqrt(2*pi*s); % functions 
if (nargout<2),return 
else 
   dfp(:,1)=(f-b0)/A; 
   dfp(:,2)=(x-a).*(f-b0)./s^2; 
   dfp(:,3)=(f-b0).*(2*x.^2-4*x*a+2*a^2-s^2)/(2*s^3); 
   dfp(:,4)=ones(size(x)); 
   dfp=dfp(:,ivar); 
end 
return 
 
function [f,dfp]=sline(x,p,ivar) 
% f=b+m*x 
x=x(:); 
b=p(1); 
m=p(2); 
f=b+m*x; 
if (nargout < 2),return 
else 
   dfp(:,2)=x; 
   dfp(:,1)=ones(size(x)); 
   dfp=dfp(:,ivar); 
end 
return 
 
function [f,dfp]=gasd(x,p,ivar) 
%f=bg+A*exp[-(x-a)^2/(2*s^2)] 
if (nargin <3),ivar=ones(length(p));end %default is vary all 
x=x(:); 
A=p(1); 
a=p(2); 
s=abs(p(3)); 
b0=p(4); 
z=(x-a)./s; 
f=b0+A*exp(-z.^2/2); % functions 
if (nargout<2),return 
else 
   dfp(:,1)=(f-b0)/A; 
   dfp(:,2)=(x-a).*(f-b0)./s^2; 
   dfp(:,3)=(f-b0).*(x-a).^2./(s^3); 
   dfp(:,4)=ones(size(x)); 
   dfp=dfp(:,ivar); 
end 
return 

function [f,dfp]=quad(x,p,ivar) 
% f=b0+b1*x + b2*x*x
x=x(:); 
bo=p(1); 
b1=p(2); 
b2=p(3);
f=bo+b1*x+b2*x.^2; 
if (nargout < 2),return 
else 
   dfp(:,1)=ones(size(x)); 
   dfp(:,2)=x;
   dfp(:,3)=x.^2;
   dfp=dfp(:,ivar); 
end 
return 

function [f,dfp]=misra1(x,p,ivar) 
%f=b1*(1-exp(-b2*x))
x=x(:); 
b1=p(1); 
b2=p(2); 
f=b1*(1-exp(-b2*x)); 
if (nargout < 2),return 
else 
   dfp(:,2)=b1*x.*exp(-b2*x); 
   dfp(:,1)=f/b1; 
   dfp=dfp(:,ivar); 
end 
return 
