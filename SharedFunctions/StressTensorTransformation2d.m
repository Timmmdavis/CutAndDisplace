function [ P11,P22,P12 ] = StressTensorTransformation2d(Pxx,Pyy,Pxy,CosAx,CosAy )
%StressTensorTransformation2d Converts 2d stress tensors into a new
%coordinate system. 

%   Copyright 2017, Tim Davis, The University of Aberdeen

%Pxx,Pyy,Pxy = input tensors 2x2 but obviously shear tensors are the same
%CosAx,CosAy = direction cosines

% %Polar 2 Cartesian stress tensors:
%[theta,R] = cart2pol(X,Y);
%CosAx=cos(theta)
%CosAy=cos((pi/2)-theta
%[ Sxx,Syy,Sxy ] = StressTensorTransformation2d(Srr,Stt,Srt,CosAx,CosAy )

% %Cartesian 2 Polar stress tensors:
%[theta,r] = cart2pol(X,-Y);
%CosAx=cos(theta)
%CosAy=cos((pi/2)-theta
%[ Srr,Stt,Srt ] = StressTensorTransformation2d(Sxx,Syy,Sxy,CosAx,CosAy);
% % http://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Polar_Coords/Polar_Coords.htm
% % -Y in cart2pol as the the martix QUAT needs to be changed, see link
% % above

%Equation 6.95 & 6.96, Pollard and Fletcher.
%Function works perfectly fine on col vectors of stress. 

%A quick low down of what is going on:
%Quat is a matrix containing directions. The first row is the first
%new coordinate direction we are transforming too. If this was a vector
%pointing along the X axis the row would be: [1,0]
%the second row is the 2nd direction. For the example above this would be
%the y axis, ie [0,1]
%Tensor is simply a matrix of the tensors in thier original input coordinates

%Unlike 3d we only need to import 1 direction vector as its easy to compute
%the other

%If the new coordinate system is a single value we repeat it. 
if numel(Pxx)~=1 %if its equal to 1 we skip as there is no point
if isscalar(CosAx) || isscalar(CosAy)
    CosAx=repmat(CosAx,size(Pxx));
    CosAy=repmat(CosAy,size(Pxx));
end    
end    
% %making sure everything is col vecs. %%Calling this is a bottleneck! 
% [ Pxx,Pyy,Pxy,CosAx,CosAy ] = RowVecToCol( Pxx,Pyy,Pxy,CosAx,CosAy );

%Superfast method:
%Assuming these are nice col vecs! 
%A var to speed stuff up:
CosAxPxx=CosAx.*Pxy;

P11q=(CosAx.*Pxx)+(-CosAy.*Pxy);
P12q=CosAxPxx+(-CosAy.*Pyy);
P21q=(CosAy.*Pxx)+CosAxPxx;
P22q=(CosAy.*Pxy)+(CosAx.*Pyy);

P11=(CosAx.*P11q)+(-CosAy.*P12q);
P12=((CosAx.*P12q)+(-CosAy.*P22q)+(CosAy.*P11q)+(CosAx.*P21q))./2;
P22=(CosAy.*P21q)+(CosAx.*P22q);



% %Preallocating blank arrays of right size 
% P11=zeros(size(Pxx));   %Prr in polar
% P22=P11;                %Ptt in polar
% P12=P11;   %Prt in polar
% 
% for i=1:numel(Pxx)
% 
% %Performing calculation for cartesian stresses
% %Eq 6.94, Pollard, arranging cosines of new directions in table
% % Quat=[CosAx(i,:),-CosAy(i,:); %dir 1, 2d vector we are transforming too
% %       CosAy(i,:),CosAx(i,:)];%dir 2, has to be 90 to 1, so is this... 
% Quat=[CosAx(i,:),-CosAy(i,:);CosAy(i,:),CosAx(i,:)];%dir 2, has to be 90 to 1, so is this... 
% 
% %Chucking the tensor together.
% % Tensor=[Pxx(i,:),Pxy(i,:);
% %         Pxy(i,:),Pyy(i,:)];
% Tensor=[Pxx(i,:),Pxy(i,:);Pxy(i,:),Pyy(i,:)];
% 
% %http://continuummechanics.org/stressxforms.html
% %Computing rotation
% CartStress=Quat*Tensor*Quat';
% 
% 
% P11(i,:)=CartStress(1,1);
% P22(i,:)=CartStress(2,2);
% P12(i,:)=CartStress(1,2);
% 
% 
% end

end
