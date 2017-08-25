function [TVA,TFA,TNA,na,TVB,TFB,TNB,nb] = stretch(VA,VB,maxarea1,maxarea2,aniso1,aniso2,E1,H1,E2,H2)
%Generates pairs of bad meshes (long triangles) to test our algorithms
%for solving the Poisson equation on the union of
%two domains.
%
%using stretch to generate bad meshes
%
%methodology: first rotate the x axis according to the degree, then
%surpress the vertices, next use function triangle to generate triangles,
%finally stretch the vertices
%
%easy to rotate in polar coordinators
%
% Inputs:
%   'VA' #VA by 2 matrix of vertices according to 'name'
%   'VB' #VB by 2 matrix of vertices according to 'name'
%   'maxArea'  a constraint on the maximum area of
%     triangles in the mesh. See the documentation
%     of the 'triangle' function.
%   'aniso1' (a,b) 1 by 2 vector where a stands for the magnitude of
%   surpression while b stands of the direction of surpression of domain 1
%   'aniso2' (a,b) 1 by 2 vector where a stands for the magnitude of
%   surpression while b stands of the direction of surpression of domain 2
%   
%
% Outputs:
%   'TVA' and 'TVB' are #TVA and #TVB by 2 matrices containing the vertices
%     of the meshes of each of the domains
%   'TFA' and 'TFB' are #TFA and #TFB by 3 matrices of indeces specifying
%     the triangular mesh
%   'TNA' and 'TNB' are #TFA and #TFB by 3
%     lists of triangle neighbour indices.
%     na (resp. nb) is a natural number such that 'TVA(1:na,:)' (resp 
%   'na' (resp. 'nb') is the number of boundary
%     vertices of 'TVA' (resp. 'TVB'). These
%     boundary vertices are listed at the
%     beginning.
    angle1 = aniso1(2)/180*pi;
    angle2 = aniso2(2)/180*pi;
    if nargin < 7
        E1 = [];
        H1 = [];
        E2 = [];
        H2 = [];
    elseif (7 < nargin) && (nargin < 9)
        if ~isempty(H1) 
            H1 = H1 * rotation(angle1);
        end
        E2 = [];
        H2 = [];
    else
        if ~isempty(H1) 
            H1 = H1 * rotation(angle1);
        end
        if ~isempty(H2)
            H2 = H2 * rotation(angle2);
        end
    end
    VA1 = VA * rotation(angle1);
    VB1 = VB * rotation(angle2);
    S1 = [aniso1(1),0;0,1];
    S2 = [1/aniso1(1),0;0,1];
    S3 = [aniso2(1),0;0,1];
    S4 = [1/aniso2(1),0;0,1];
    VA1 = VA1 * S1;
    VB1 = VB1 * S3;
    if isempty(E1)
        [TVA1,TFA,TNA,na]=triangle_bd_int(VA1,'NoBoundarySteiners','MaxArea',aniso1(1)* maxarea1);
    else
        [TVA1,TFA,TNA,na]=triangle_bd_int(VA1,E1,H1,'NoBoundarySteiners','MaxArea',aniso1(1)* maxarea1);
    end
    if isempty(E2)
        [TVB1,TFB,TNB,nb]=triangle_bd_int(VB1,'NoBoundarySteiners','MaxArea',aniso2(1)* maxarea2);
    else
        [TVB1,TFB,TNB,nb]=triangle_bd_int(VB1,E2,H2,'NoBoundarySteiners','MaxArea',aniso2(1)*maxarea2);
    end
    TVA1 = TVA1 * S2;
    TVB1 = TVB1 * S4;
    TVA = TVA1 * rotation(-angle1);
    TVB = TVB1 * rotation(-angle2);
end

