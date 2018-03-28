function [ FaceNormalVector ] = CalculateTriangleNormal( Pa,Pb,Pc )
%CalculateTriangleNormal Calculates the triangle normal from 3 input
%points (edge corners of triangle). 
%Pa=(X1,Y1,Z1)
%Pb=(X2,Y2,Z2)
%Pc=(X3,Y3,Z3)

%Now calculating the normal orientation
%New vectors (see calculating normals online)
U=Pb-Pa;
V=Pc-Pa;
%Cross product of the vectors
Nx = (U(2)*V(3)) - (U(3)*V(2));
Ny = (U(3)*V(1)) - (U(1)*V(3));
Nz = (U(1)*V(2)) - (U(2)*V(1));
%Vector Magnitude
aMag=sqrt((Nx * Nx) + (Ny * Ny) + (Nz * Nz));
%Norm values
Ax=Nx/aMag;
Ay=Ny/aMag;
Az=Nz/aMag;
FaceNormalVector=[Ax,Ay,Az];

end

