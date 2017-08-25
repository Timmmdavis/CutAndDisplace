function [ P11,P22,P33,P12,P13,P23 ] = StressTensorTransformation3d(Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,C1,C2,C3)
%StressTensorTransformation3d Converts 3d stress tensors into a new
%coordinate system. 

%Pxx,Pyy,Pxy etc are the tensors in the current coordinate system
%C1 C2 C3 are 3*n vector the direction cosines of the new
%coordinate system
%P11,P22,P33,P12,P13,P23 = output tensors that are in the new coords

%   Copyright 2017, Tim Davis, The University of Aberdeen

%Does Equations 6.88 - 6.93 of Pollard and Fletcher in one step, 
%Function works perfectly fine on col vectors of stress. 

%A quick low down:
%Quat is a matrix containing directions (quaternary). The first row is the first
%new coordinate direction we are transforming too. If this was a vector
%pointing along the X axis the row would be: [1,0,0] or [A11,A12,A13] as it
%only points along X and is 90 to the other two directions
%the second row is the 2nd direction. For the example above this would be
%the y axis, ie [0,1,0]
%etc
%Tensor is simply a matrix of the tensors in their original input coordinates

%If the new coordiante system is a single value we repeat it. 
if isscalar(C1(:,1)) || isscalar(C2(:,1)) || isscalar(C3(:,1)) %checks if its just 1 row
    C1=repmat(C1,size(Pxx(:,1)));
    C2=repmat(C2,size(Pxx(:,1)));
    C3=repmat(C3,size(Pxx(:,1)));
end  

A11=C1(:,1);A12=C1(:,2);A13=C1(:,3);
A21=C2(:,1);A22=C2(:,2);A23=C2(:,3);
A31=C3(:,1);A32=C3(:,2);A33=C3(:,3);
    
%making sure everything is col vecs
[ Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,A11,A12,A13,A21,A22,A23,A31,A32,A33 ] = RowVecToCol( Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,A11,A12,A13,A21,A22,A23,A31,A32,A33 );

for i = 1:numel(Pxx)
%checking all vectors are 90 to each other
%this should give a result of 0 if there is an angle of 90 between the vectors 
DotPro12=A11(i).*A12(i)+A21(i).*A22(i)+A31(i).*A32(i);
DotPro13=A11(i).*A13(i)+A21(i).*A23(i)+A31(i).*A33(i);
DotPro23=A12(i).*A13(i)+A22(i).*A23(i)+A32(i).*A33(i);
    if any(round(DotPro12,14)) || any(round(DotPro13,14)) || any(round(DotPro23,14)) %rounding so no precision issues close to eps
        error('Direction cosines are not at 90 degrees to each other')
    end    
end

%Preallocating blank arrays of right size 
P11=zeros(size(Pxx));   %in cart xx
P22=P11;                %in cart yy
P33=P11;                %in cart zz
P12=P11;                %in cart xy
P13=P11;                %in cart xz
P23=P11;                %in cart yx


for i=1:numel(Pxx)

% Performing calculation for cartesian stresses
%Eq 2.23, Pollard, arranging cosines in table
Quat=[A11(i,:),A12(i,:),A13(i,:); %in cart dir cosines of the x axis
      A21(i,:),A22(i,:),A23(i,:); %in cart y axis
      A31(i,:),A32(i,:),A33(i,:)];%in cart z axis

%Chucking the tensor together, if you have somehow managed to interpret/compute shear components of
%stress from field observations you could put them in here too
Tensor=[Pxx(i,:),Pxy(i,:),Pxz(i,:);
        Pxy(i,:),Pyy(i,:),Pyz(i,:);
        Pxz(i,:),Pyz(i,:),Pzz(i,:)];

%http://continuummechanics.org/stressxforms.html
%Computing rotation
CartStress=Quat*Tensor*Quat';

% Extracting variables
P11(i,:)=CartStress(1,1);
P22(i,:)=CartStress(2,2);
P33(i,:)=CartStress(3,3);
P12(i,:)=CartStress(1,2);
P13(i,:)=CartStress(1,3);
P23(i,:)=CartStress(2,3);
end

end

