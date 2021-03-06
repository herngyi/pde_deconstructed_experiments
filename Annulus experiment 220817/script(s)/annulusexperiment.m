%   EXPERIMENT II - (DATE OF UPLOAD)
%We try out all of our methods (including remeshing) on two concentric intersecting annulus for
%different resolutions (represented as average edge length), and compare
%them against the known analytic solution.
%
g=@(V) (sqrt((V(:,1).^2)+(V(:,2).^2)))<1.5;
innerbc=1;
outerbc=0;
A=[1,log(1);1,log(2)];
b=[1;0];
C=inv(A)*b;
C1=C(1);
C2=C(2);
gt =@(V) C1+C2*log(sqrt((V(:,1).^2)+(V(:,2).^2)));
errordelta=[];
errornaive=[];
errorugt=[];
%errorneumann=[];
h=[];
hgt=[];
for s=5:11
    utrue=[];
    [VA,FA,NA]=annulus(2^s,2,'R',1.4);
    [VB,FB,NB]=annulus(2^(s-1),1.6,'R',1); 
    [V,F]=annulus(2^s,2,'R',1);
    v=outline(F);
    v=unique(v(:));
    ZZ=overlap_poisson({VA,VB},{FA,FB},g,@(V) zeros(size(V,1),1),'Method','dirichlet');
    udelta=[ZZ{1};ZZ{2}];
    ZZ=overlap_poisson({VA,VB},{FA,FB},g,@(V) zeros(size(V,1),1),'Method','naive');
    unaive=[ZZ{1};ZZ{2}];
    Q=cotmatrix(V,F);
    ugt=min_quad_with_fixed(Q,sparse(size(Q,1),1),v,g(V(v,:)));
    eugt=[];
    %[uAneumann,uBneumann]=solve_intersecting(VA,FA,NA,va,VB,FB,NB,vb,g,'neumann');
    %uneumann=[uAneumann;uBneumann];
     utrue = gt([VA;VB]);
    utruegt = gt(V);
    errordelta=[errordelta,max(abs(utrue-udelta))];
    errornaive=[errornaive,max(abs(utrue-unaive))];
    errorugt=[errorugt,max(abs(ugt-utruegt))];
    %errorneumann=[errorneumann,max(abs(utrue-uneumann))];
    disp(length(errordelta))
    h=[h,avgedge(VB,FB)];
    hgt =[hgt,avgedge(V,F)];
    H= [h',h',hgt'];
    E = [errordelta',errornaive',errorugt'];
loglog(H,E,'LineWidth',3)
    axis equal
   % legend('DSC','OSC','GT')
    title('Convergence for annulus test')
    xlabel('h')
    ylabel('L_\infty error')
    drawnow
    saveas(gcf,'annulus2dexp','epsc')
end
