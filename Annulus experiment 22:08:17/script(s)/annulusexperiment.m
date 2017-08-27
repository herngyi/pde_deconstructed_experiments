%   EXPERIMENT II - (DATE OF UPLOAD)
%We try out all of our methods (including remeshing) on two concentric intersecting annulus for
%different resolutions (represented as average edge length), and compare
%them against the known analytic solution.
%
g=@(x,y) (sqrt((x.^2)+(y.^2)))<1.5
innerbc=1;
outerbc=0;
A=[1,log(1);1,log(2)];
b=[1;0];
C=inv(A)*b;
C1=C(1);
C2=C(2);
errordelta=[];
errornaive=[];
errorugt=[];
%errorneumann=[];
h=[];
for s=3:10
    utrue=[];
    [VA,FA,NA]=annulus_neigh(2^s,2,'R',1.4);
    [VA,FA,va]=bd_loops_first(VA,FA);
    [VB,FB,NB]=annulus_neigh(2^s,1.6,'R',1);
    [VB,FB,vb]=bd_loops_first(VB,FB);   
    TH=[0,0];
    [V,F,N,v] = remesh_union(VA,FA,va,VB,FB,vb,TH);
    [uAdelta,uBdelta]=solve_intersecting(VA,FA,NA,va,VB,FB,NB,vb,g,'delta');
    udelta=[uAdelta;uBdelta];
    [uAnaive,uBnaive]=solve_intersecting(VA,FA,NA,va,VB,FB,NB,vb,g,'naive');
    unaive=[uAnaive;uBnaive];
    Q=cotmatrix(V,F);
    ugt=min_quad_with_fixed(Q,sparse(size(Q,1),1),1:v(end),g(V(1:v(end),1),V(1:v(end),2)));
    eugt=[];
    %[uAneumann,uBneumann]=solve_intersecting(VA,FA,NA,va,VB,FB,NB,vb,g,'neumann');
    %uneumann=[uAneumann;uBneumann];
    for i=1:size(VA,1)
        utrue=[utrue;C1+C2*log(sqrt((VA(i,1)^2)+(VA(i,2)^2)))];
        eugt=[eugt;abs(ugt(i)-C1-C2*log(sqrt((V(i,1)^2)+(V(i,2)^2))))];
    end
    for i=1:size(VB,1)
        utrue=[utrue;C1+C2*log(sqrt((VB(i,1)^2)+(VB(i,2)^2)))];
        eugt=[eugt;abs(ugt(i+size(VA,1))-C1-C2*log(sqrt((V(i+size(VA,1),1)^2)+(V(i+size(VA,1),2)^2))))];
    end
    errordelta=[errordelta,max(abs(utrue-udelta))];
    errornaive=[errornaive,max(abs(utrue-unaive))];
    errorugt=[errorugt,max(eugt)];
    %errorneumann=[errorneumann,max(abs(utrue-uneumann))];
    disp(length(errordelta))
    h=[h,avgedge(VA,FA)];
end
plot(log(h),log(errordelta),log(h),log(errornaive),log(h),log(errorugt))
axis equal
legend('DSC','OSC','GT')
title('Convergence for annulus test')
xlabel('log avg edge length')
ylabel('log linf error')
