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
    %annulus with inner radius 1.4 and outer radius 2
    region1 = circle(0,0,2);
    region2 = circle(0,0,1.4);
    %annulus with inner radius 1 and out radius 1.6
    region3 = circle(0,0,1.6);
    region4 = circle(0,0,1);
    A = [region1;region2];
    B = [region3;region4];
    a(1) = size(region1,1);
    a(2) = size(region2,1);
    b(1) = size(region3,1);
    b(2) = size(region4,1); 
    H1 = [0,0];
    E1 = [(a(1)+1):(a(1)+a(2)),1:a(1);(a(1)+2):(a(1)+a(2)),(a(1)+1),2:a(1),1]';
    H2 = [0,0];
    E2 = [(b(1)+1):b(1)+b(2),1:b(1);(b(1)+2):(b(1)+a(2)),(b(1)+1),2:b(1),1]';
    [VA,FA,NA,va,VB,FB,NB,vb] = stretch(A,B,maxarea,maxarea,aniso,aniso,E1,H1,E2,H2);
    TH=[0,0];
    [V,F,N,cycLocs] = remesh_union(VA,FA,va,VB,FB,vb,TH);
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
