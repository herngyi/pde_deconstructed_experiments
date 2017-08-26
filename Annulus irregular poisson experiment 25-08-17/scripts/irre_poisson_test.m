function [errordelta,errornaive,errorugt,h] = irre_poisson_test(maxarea,aniso)
g=@(x,y) x.*0;
H=@(x,y) [-4.*x,x.*0]
R0=1
R1=2
poisson_exact = @(r) -r.^2+(R0.^2-R1.^2)./(log(R0)-log(R1)).*log(r)+(log(R0).*R1.^2-R0.^2.*log(R1))./(log(R0)-log(R1));
errordelta=[];
errornaive=[];
errorugt=[];
%errorneumann=[];
h=[];
%generate the bad meshes
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
    [TV,TF,N,cycLocs] = remesh_union(VA,FA,va,VB,FB,vb,TH);
    utrue=[];
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
  
    for i=1:size(VA,1)
        utrue=[utrue;poisson_exact(sqrt((VA(i,1)^2)+(VA(i,2)^2)))];
        eugt=[eugt;abs(ugt(i)-poisson_exact(sqrt((TV(i,1)^2)+(TV(i,2)^2))))];
    end
    for i=1:size(VB,1)
        utrue=[utrue;poisson_exact(sqrt((VB(i,1)^2)+(VB(i,2)^2)))];
        eugt=[eugt;abs(ugt(i+size(VA,1))-poisson_exact(sqrt((TV(i+size(VA,1),1)^2)+(TV(i+size(VA,1),2)^2))))];
    end
    errordelta= max(abs(utrue-udelta));
    errornaive= max(abs(utrue-unaive));
    errorugt= max(eugt);
    h=[h,avgedge(VA,FA)];