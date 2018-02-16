%   EXPERIMENT II - (DATE OF UPLOAD)
%We try out all of our methods (including remeshing) on two concentric intersecting annulus for
%different resolutions (represented as average edge length), and compare
%them against the known analytic solution of the poisson equation with
%constant right-hand side.
%
g=@(x,y) x.*0;
H=@(x,y) [x,x.*0]
R0=1
R1=2
poisson_exact = @(r) (-1/4).*(-r.^2+(R0.^2-R1.^2)./(log(R0)-log(R1)).*log(r)+(log(R0).*R1.^2-R0.^2.*log(R1))./(log(R0)-log(R1)));
errordelta=[];
errornaive=[];
errorugt=[];
%errorneumann=[];
h=[];
for s=3:12
    utrue=[];
    [VA,FA,NA]=annulus_neigh(2^s,2,'R',1.2);
    [VA,FA,va]=bd_loops_first(VA,FA);
    [VB,FB,NB]=annulus_neigh(2^s,1.8,'R',1);
    [VB,FB,vb]=bd_loops_first(VB,FB);   
    TH=[0,0];
    [TV,TF,N,cycLocs] = remesh_union(VA,FA,va,VB,FB,vb,TH);
    G=grad(TV,TF);
    dbl=doublearea(TV,TF);
    dbl=repdiag(diag(sparse(dbl)/2),size(TV,2));
    X=H((TV(TF(:,1),1)+TV(TF(:,2),1)+TV(TF(:,3),1))./3,(TV(TF(:,1),2)+TV(TF(:,2),2)+TV(TF(:,3),2))./3);
    X=[X(:,1);X(:,2)];
    Q=G'*dbl*G;
    B=(-X'*dbl*G-(G'*dbl*X)')';
    ugt=min_quad_with_fixed(Q,B,[1:cycLocs(end)]',g(1:cycLocs(end))');
    [uAdelta,uBdelta]=solve_intersecting_poisson(VA,FA,NA,va,VB,FB,NB,vb,g,H,'delta');
    udelta=[uAdelta;uBdelta];
    [uAnaive,uBnaive]=solve_intersecting_poisson(VA,FA,NA,va,VB,FB,NB,vb,g,H,'naive');
    unaive=[uAnaive;uBnaive];
    
    
    eugt=[];
    %[uAneumann,uBneumann]=solve_intersecting(VA,FA,NA,va,VB,FB,NB,vb,g,'neumann');
    %uneumann=[uAneumann;uBneumann];
    for i=1:size(VA,1)
        utrue=[utrue;poisson_exact(sqrt((VA(i,1)^2)+(VA(i,2)^2)))];
        eugt=[eugt;abs(ugt(i)-poisson_exact(sqrt((TV(i,1)^2)+(TV(i,2)^2))))];
    end
    for i=1:size(VB,1)
        utrue=[utrue;poisson_exact(sqrt((VB(i,1)^2)+(VB(i,2)^2)))];
        eugt=[eugt;abs(ugt(i+size(VA,1))-poisson_exact(sqrt((TV(i+size(VA,1),1)^2)+(TV(i+size(VA,1),2)^2))))];
    end
    errordelta=[errordelta,max(abs(utrue-udelta))];
    errornaive=[errornaive,max(abs(utrue-unaive))];
    errorugt=[errorugt,max(eugt)];
    %errorneumann=[errorneumann,max(abs(utrue-uneumann))];
    disp(length(errordelta))
    h=[h,avgedge(VA,FA)];
end
plot(log(h),log(errornaive),log(h),log(errordelta),log(h),log(errorugt))
axis equal
legend('OSC','DSC','GT')
title('Convergence for annulus test')
xlabel('log avg edge length')
ylabel('log linf error')