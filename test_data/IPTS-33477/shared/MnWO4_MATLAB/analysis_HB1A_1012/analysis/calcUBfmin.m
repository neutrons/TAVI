function [estimates, model] = calcUBfmin(data, pa, ia)

    idx=find(ia);
    start_point=pa(idx);
    model = @ubmatrixmin;
    estimates = fminsearch(model, start_point,optimset('MaxFunEvals',3e4,'MaxIter',5000,'TolX',1e-5,'Display','iter'));

    function [sse,fval,sample, UBmatrix] = ubmatrixmin(p)
    %initialize the parameters to be used in matrix calculation.
    idx=find(ia);
    newp=zeros(11,1);
    for i=1:length(idx)
	newp(idx(i))=p(i);
    end
    idx=find(ia==0);
    for i=1:length(idx)
	newp(idx(i))=pa(idx(i));
    end
    p=newp(:);

    lambda=1.536;
    %initialize the UBmatrix
    UB11=p(1); UB12=p(2); UB13=p(3);
    UB21=p(4); UB22=p(5); UB23=p(6);
    UB31=p(7); UB32=p(8); UB33=p(9);
    theta20=p(10); chi0=p(11);

    UBmatrix=[UB11,UB12,UB13;UB21,UB22,UB23;UB31,UB32,UB33];
    UBtilt=UBmatrix';
    InvG=UBtilt*UBmatrix;
    G=inv(InvG);

    B=[];
    qi=[];
    sample.a=sqrt(G(1,1));
    sample.b=sqrt(G(2,2));
    sample.c=sqrt(G(3,3));
    sample.alpha=acos(G(2,3)/sample.b/sample.c);
    sample.beta=acos(G(1,3)/sample.a/sample.c);
    sample.gamma=acos(G(1,2)/sample.a/sample.b);
    [B,V,Vstar,latticestar]=Bmatrix(sample);

    h=data(:,1);
    k=data(:,2);
    l=data(:,3);
    theta2=(data(:,4)-theta20)*pi/180;
    omega=(data(:,5)-data(:,4)/2)*pi/180;
    chi=(data(:,6)-chi0)*pi/180;
    phi=data(:,7)*pi/180;

    Q=[h(:) k(:) l(:)]';
    TT=UBmatrix*Q;
     
    % Calculate the matrix of u_phi
    U(1,:)=cos(omega).*cos(chi).*cos(phi)-sin(omega).*sin(phi);
    U(2,:)=cos(omega).*cos(chi).*sin(phi)+sin(omega).*cos(phi);
    U(3,:)=cos(omega).*sin(chi);

    sse=0;
    qi=2*sin(theta2(:)/2)'./lambda;
    err=TT-(repmat(qi,3,1).*U);
    errsq=err.*err;
    sse=sum(errsq(:));
    hphi=repmat(qi,3,1).*U;
    fval=inv(UBmatrix)*hphi;
end
end
