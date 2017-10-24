function [S1,S2,S3,S1dir,S2dir,S3dir] = EigCalc3dSpeed(Sxx,Syy,Szz,Sxy,Sxz,Syz)
% Calculates the 3-D principal stress/strain magnitudes and directions from input tensors
% Input parameters
% Stress or strain tensors 9x9 but Sxy=Syx so we only bring in 6. 
% Output parameter (tension positive convention for inputs)
% S1 = most tensile principal stress
% S3 = most comp principal stress
% S2 = intermediate
% S1dir = XYZ components of the direction of the principal stress S1
% S2dir & S3dir = see above
% Examples:
%Calclating 3d stress eigenValues
%[S1,S2,S3,S1dir,S2dir,S3dir] = EigCalc3d(Sxx(:),Syy(:),Szz(:),Sxy(:),Sxz(:),Syz(:));
%Calclating 2d strain eigenValues
%[E1,E2,E3,E1dir,E2dir,E3dir] = EigCalc3d(Exx(:),Eyy(:),Ezz(:),Exy(:),Exz(:),Eyz(:));
%   Copyright 2017, Tim Davis, The University of Potsdam

%Preallocating array
n = numel(Sxx);
tensor=zeros (3,3,n);

%Filling this array with 3x3 tensors, each 3rd dim is each point
%Accumulated using indexing
%[Sxx,Sxy,Sxz]
%[Syx,Syy,Syz]
%[Szx,Szy,Szz]
tensor(1,1,1:1:end) = Sxx(1:1:end,:);
tensor(1,2,1:1:end) = Sxy(1:1:end,:);
tensor(1,3,1:1:end) = Sxz(1:1:end,:);
tensor(2,1,1:1:end) = Sxy(1:1:end,:);
tensor(2,2,1:1:end) = Syy(1:1:end,:);
tensor(2,3,1:1:end) = Syz(1:1:end,:);
tensor(3,1,1:1:end) = Sxz(1:1:end,:);
tensor(3,2,1:1:end) = Syz(1:1:end,:);
tensor(3,3,1:1:end) = Szz(1:1:end,:);

%Eig can't handle nan's so we turn these to 0's and put the calculated s1s2s3 to nans after 
NanFlag = isnan(tensor);
tensor(NanFlag)=0;

%Do the calculation
D = eig3(tensor); %see base of script and functions it calls
%Sort results
D = sortrows(D);
%Transpose
D=D';
%Extract
S3 = D(:,1);
S2 = D(:,2);
S1 = D(:,3);

%Create vars for loop (preallocate)
S1dir=zeros(3,n);
S2dir=zeros(3,n);
%Create Identity matrix
Ident = repmat(eye(3),1,1,n);

for i=1:n
    %Equation we want to solve
    %[A-(Lambda * Identity)]v=0
    
    %Grab the current bit in the loop
    I=Ident(:,:,i);
    A=tensor(:,:,i);
    
    %[Lambda * Identity] in eigenvector equation
    One=S1(i)*I;
    Two=S2(i)*I;
    
    %[A-(Lambda * Identity)]
    vect1=A-One;
    vect2=A-Two;
    
    %[A-(Lambda * Identity)]v=0. Finding vector [v]
    %https://de.mathworks.com/help/matlab/ref/decomposition.html
    [Q1,~] = qr(vect1,0);
    [Q2,~] = qr(vect2,0);

    %Extracting results
    S1dir(:,i)=(Q1(:,3));
    S2dir(:,i)=(Q2(:,3));
    
end

S1dir=S1dir';
S2dir=S2dir';

%Doing cross on these big arrays to get S3 direction, we know its
%perpendicular to S2 and S1. 
S3dir=cross2(S1dir,S2dir);

end

function [D]=eig3(A)
% function D = eig3(A)
% 
% Copyright (c) 2010, Bruno Luong
% All rights reserved.
%
%
% Compute in one shot the eigen-values of multiples (3 x 3) matrices
%
% INPUT:
%   A: (3 x 3 x n) array
% OUTPUT:
%   D: (3 x n). EIG3 returns in D(:,k) three eigen-values of A(:,:,k)
%
% See also: CardanRoots, eig2, eig
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History:
%     Original 20-May-2010

if size(A,1) ~= 3 || size(A,2) ~= 3
    error('A must be [3x3xn] array');
end

A = reshape(A, 9, []).';

P3 = 1;
% Trace
P2 = -(A(:,1)+A(:,5)+A(:,9));
% Principal minors
M11 = A(:,5).*A(:,9) - A(:,8).*A(:,6);
M22 = A(:,9).*A(:,1) - A(:,3).*A(:,7);
M33 = A(:,1).*A(:,5) - A(:,4).*A(:,2);
P1 = (M11 + M22 + M33);
% Determinant
P0 = - A(:,1).*M11 ...
     + A(:,4).*(A(:,2).*A(:,9)-A(:,8).*A(:,3)) ...
     - A(:,7).*(A(:,2).*A(:,6)-A(:,5).*A(:,3));

D = CardanRoots(P3, P2, P1, P0).';

end