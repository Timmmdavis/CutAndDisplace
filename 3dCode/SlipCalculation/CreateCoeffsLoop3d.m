function [infmatrix,Dinfmatrix]=CreateCoeffsLoop3d(infmatrix,Dinfmatrix,...
    NUM,X,Y,Z,P1,P2,P3,Dss,Dds,Dn,mu,lambda,nu,halfspace,FD)
% CreateCoeffsLoop3d: Loop that calls the TDE functions of Mehdi 
%               Nikkhoo and fills a large column mat with influence
%               coefficients.
%
% usage #1:
% [infmatrix,DisplacementXYZ]=CreateCoeffsLoop3d(infmatrix,DisplacementXYZ,...
%     NUM,X,Y,Z,P1,P2,P3,Dss,Dds,Dn,mu,lambda,nu,halfspace,FD)
%
% Arguments: (input)
%  infmatrix  - 6*n 'vector' of 0's to be filled with strain coefficients.
%
% Dinfmatrix  - 3*n 'vector' of 0's to be filled with displacement coefficients.
%
%       NUM   - The number of elements in the calculation. 
%
%    X,Y,Z    - The element midpoints (triangle). (col vec, XYZ)
%
% P1,P2,P3    - n*3 Column vectors where each 'P' represents the
%               different corner points of one of the triangles (XYZ).
%
% Dss,Dds,Dn  - One will be one and the others 0, depending if you are creating
%              a Dss,Dds or Dn inf matrix.
%
%       mu    - Shear modulus.
%
%     lambda  - Lame's constant.
%
%       nu    - The Poisson's ratio.
%
%  halfspace  - Defines if we work out the coefficientsin a half or whole
%              space
%
%
% FD          - Flag telling the user if any elements are going to be fixed
%              (if this is the case displacement influence matricies are
%              filled).
%
% Arguments: (output)
%   infmatrix  - Strains Exx,Eyy,Ezz,Exy,Exz,Eyz as a 'vector' where
%               every NUM rows in the vector is the influence of a new
%               element. 
%
%   Dinfmatrix  - Displacement Dx,Dy,Dz as a 'vector' where
%               every NUM rows in the vector is the influence of a new
%               element.   
%
% Example usage:
%
% [infmatrix,DisplacementXYZ]=CreateCoeffsLoop3d(infmatrix,DisplacementXYZ,...
%     NUM,X,Y,Z,P1,P2,P3,Dss,Dds,Dn,mu,lambda,nu,halfspace,FD)
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%2nd size 
NUM2=size(P1,1);

%Strings for progress
if Dss==1
    StringHS='1/3 CalculatingDssInfMatrixHS';
    StringFS='1/3 CalculatingDssInfMatrixFS';
elseif Dds==1
    StringHS='2/3 CalculatingDdsInfMatrixHS';
    StringFS='2/3 CalculatingDdsInfMatrixFS';
else %Dn==1
    StringHS='3/3 CalculatingDnInfMatrixHS';
    StringFS='3/3 CalculatingDnInfMatrixFS';
end

%Starting progress bar function
if halfspace==1
    %Getting progress bar string
    progressbar(StringHS)  
else
    %Getting progress bar string
    progressbar(StringFS)          
end       


for i=1:size(P1,1)
    %Creating size and space that will be filled in the influence
    %matrix in each loop
    first = (i-1)*NUM+1;
    last = i*NUM;
    if halfspace==1
        infmatrix(first:last,:) = TDstrainHS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Dss,Dds,Dn,mu,lambda);
        if FD==1
            Dinfmatrix(first:last,:) = TDdispHS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Dss,Dds,Dn,nu);
        end
    else %Fullspace
        infmatrix(first:last,:) = TDstrainFS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Dss,Dds,Dn,mu,lambda);
        if FD==1
            Dinfmatrix(first:last,:) = TDdispFS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Dss,Dds,Dn,nu); 
        end
    end
    %Updating progress bar
    progressbar(i/NUM2) % Update figure 
end



end