function [errordelta,errornaive,errorugt,h] = irre_converge_test(maxarea,aniso)
g=@(x,y) (sqrt((x.^2)+(y.^2)))<1.5;
gg=@(V) (sqrt((V(:,1).^2)+(V(:,2).^2)))<1.5;
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
h=[];

%generate the bad meshes
    utrue=[];
    %annulus with inner radius 1.4 and outer radius 2
    region1 = circle(0,0,2);
    region2 = circle(0,0,1.4);
    %annulus with inner radius 1 and out radius 1.6
    region3 = circle(0,0,1.6);
    region4 = circle(0,0,1);
    VA = [region1;region2];
    VB = [region3;region4];
    va(1) = size(region1,1);
    va(2) = size(region2,1);
    vb(1) = size(region3,1);
    vb(2) = size(region4,1); 
    H1 = [0,0];
    E1 = [(va(1)+1):(va(1)+va(2)),1:va(1);(va(1)+2):(va(1)+va(2)),(va(1)+1),2:va(1),1]';
    H2 = [0,0];
    E2 = [(vb(1)+1):vb(1)+vb(2),1:vb(1);(vb(1)+2):(vb(1)+va(2)),(vb(1)+1),2:vb(1),1]';
    [VA,FA,NA,na,VB,FB,NB,nb] = stretch(VA,VB,maxarea,maxarea,aniso,aniso,E1,H1,E2,H2);
    TH=[0,0];
    [V,F,N,v] = remesh_union(VA,FA,na,VB,FB,nb,TH);
    % using dirichlet constraint on bad meshes
    ZZ=overlap_poisson({VA,VB},{FA,FB},gg,@(V) zeros(size(V,1),1),'Method','dirichlet');
    udeltabad=[ZZ{1};ZZ{2}];
    ZZ=overlap_poisson({VA,VB},{FA,FB},gg,@(V) zeros(size(V,1),1),'Method','naive');
    unaivebad=[ZZ{1};ZZ{2}];
    % using corefined good meshes - ground truth
    Q=cotmatrix(V,F);
    ugtbad=min_quad_with_fixed(Q,sparse(size(Q,1),1),1:v(end),g(V(1:v(end),1),V(1:v(end),2)));
    utruebad = [];
    eugtbad=[];
  
    for i=1:size(VA,1)
        utruebad=[utruebad;C1+C2*log(sqrt((VA(i,1)^2)+(VA(i,2)^2)))];
        eugtbad=[eugtbad;abs(ugtbad(i)-C1-C2*log(sqrt((V(i,1)^2)+(V(i,2)^2))))];
    end
    for i=1:size(VB,1)
        utruebad=[utruebad;C1+C2*log(sqrt((VB(i,1)^2)+(VB(i,2)^2)))];
        eugtbad=[eugtbad;abs(ugtbad(i+size(VA,1))-C1-C2*log(sqrt((V(i+size(VA,1),1)^2)+(V(i+size(VA,1),2)^2))))];
    end
    errordelta= max(abs(utruebad-udeltabad));
    errornaive= max(abs(utruebad-unaivebad));
    errorugt= max(eugtbad);
    h=[h,avgedge(VA,FA)];
end
