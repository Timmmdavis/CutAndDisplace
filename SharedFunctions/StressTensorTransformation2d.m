function [ P11,P22,P12 ] = StressTensorTransformation2d(Pxx,Pyy,Pxy,CosAx,CosAy )
% StressTensorTransformation2d: Converts 2d stress/strain tensors into a new
%                   coordinate system. Can also be used to convert tensors
%                   between polar and Cartesian coordinates.
%                   %Equation 6.95 & 6.96, Pollard and Fletcher.
%
%                   Good reference: http://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Polar_Coords/Polar_Coords.htm
%
% usage #1: Just changing the coordinate axis orientation. 
% [ P11,P22,P12 ] = StressTensorTransformation2d(Pxx,Pyy,Pxy,CosAx,CosAy )
%
% usage #2: Convert Polar to Cart tensors:
% [theta,R] = cart2pol(X,Y);
% CosAx=cos(theta)
% CosAy=sin(theta)
% [ Sxx,Syy,Sxy ] = StressTensorTransformation2d(Srr,Stt,Srt,CosAx,CosAy )
%
% usage #3: Convert Cart to Polar tensors:
% [theta,r] = cart2pol(X,-Y);
% CosAx=cos(theta)
% CosAy=sin(theta)
% [ Srr,Stt,Srt ] = StressTensorTransformation2d(Sxx,Syy,Sxy,CosAx,CosAy);
%
% Arguments: (input)
% Pxx,Pyy,Pxy       - 2D stress tensors (or strain). (Col vectors or mat)
%
% CosAx,CosAy       - Direction cosines of the new coordiantes for each of
%                     the tensors. (Col vectors or mat)
%
% Arguments: (output)
% P11,P22,P12       - The new tensors. (Col vectors or mat)
% 
%
% Example usage:
% 
%  [X,Y]=meshgrid(-2:0.1:2,-2:0.1:2);
%  a = 1;
%  P = 1; 
%  Sxx = -1; 
%  Syy = -1;
%  [Sxx,Syy,Sxy,Srr,Stt,Srt]=Kirsch1898_PressurisedHole(X,Y,a,P,Sxx,Syy);
%  % Radial coordiantes for the input observation point locations
%  [theta,r] = cart2pol(X,-Y);
%  CosAx=cos(theta);
%  CosAy=sin(theta);
%  %Now converting the Cart tensors to radial and comparing those to the
%  %original Kirsch function. 
%  [ SrrRotTest,SttRotTest,SrtRotTest ]...
%  = StressTensorTransformation2d(Sxx,Syy,Sxy,CosAx,CosAy);
%  SrrRotTest=reshape(SrrRotTest,41,41);
%  SttRotTest=reshape(SttRotTest,41,41);
%  SrtRotTest=reshape(SrtRotTest,41,41);
%  DrawContourFPlots2d( X,Y,[],Srr,SrrRotTest,Stt,SttRotTest,Srt,SrtRotTest );
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen



%Unlike 3d we only need to import 1 direction vector as its easy to compute
%the other

%If the new coordinate system is a single value we repeat it. 
if numel(Pxx)~=1 %if its equal to 1 we skip as there is no point
    if isscalar(CosAx) || isscalar(CosAy)
        CosAx=repmat(CosAx,size(Pxx));
        CosAy=repmat(CosAy,size(Pxx));
    end    
end    


%Superfast method:
%First make a variable to speed stuff up:
CosAxPxx=CosAx.*Pxy;

%Just doing matrix multiplication (shown below) but with indexing. 
P11q=(CosAx.*Pxx)+(-CosAy.*Pxy);
P12q=CosAxPxx+(-CosAy.*Pyy);
P21q=(CosAy.*Pxx)+CosAxPxx;
P22q=(CosAy.*Pxy)+(CosAx.*Pyy);

P11=(CosAx.*P11q)+(-CosAy.*P12q);
P12=((CosAx.*P12q)+(-CosAy.*P22q)+(CosAy.*P11q)+(CosAx.*P21q))./2;
P22=(CosAy.*P21q)+(CosAx.*P22q);


% %Old slow version:
% %A quick low down of what is going on:
% %Quat is a matrix containing directions. The first row is the first
% %new coordinate direction we are transforming too. If this was a vector
% %pointing along the X axis the row would be: [1,0]
% %the second row is the 2nd direction. For the example above this would be
% %the y axis, ie [0,1]
% %Tensor is simply a matrix of the tensors in thier original input coordinates
%
% %Preallocating blank arrays of right size 
% P11=zeros(size(Pxx));   %Prr in polar
% P22=P11;                %Ptt in polar
% P12=P11;   %Prt in polar
% 
% for i=1:numel(Pxx)
% 
%     %Performing calculation for cartesian stresses
%     %Eq 6.94, Pollard, arranging cosines of new directions in table
%     Quat=[[CosAx(i,:),-CosAy(i,:)]  %dir 1, 2d vector we are transforming too
%           [CosAy(i,:), CosAx(i,:)]];%dir 2, has to be 90 to 1, so is this... 
% 
%     %Chucking the tensor together.
%     Tensor=[[Pxx(i,:),Pxy(i,:)]
%             [Pxy(i,:),Pyy(i,:)]];
% 
%     %http://continuummechanics.org/stressxforms.html
%     %Computing rotation
%     CartStress=Quat*Tensor*Quat';
% 
%     P11(i,:)=CartStress(1,1);
%     P22(i,:)=CartStress(2,2);
%     P12(i,:)=CartStress(1,2);
% 
% end

end
