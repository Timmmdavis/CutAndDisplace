function [S1,S2,S3,S1dir,S2dir,S3dir] = EigCalc3d(Sxx,Syy,Szz,Sxy,Sxz,Syz)
% EigCalc3d: Calculates the 3D principal stress/strain magnitudes and 
%                   directions from input tensors.
%   
% usage #1: For stress:
% [S1,S2,S3,S1dir,S2dir,S3dir] = EigCalc3d(Sxx,Syy,Szz,Sxy,Sxz,Syz)
%
% usage #2: For strain:
% [E1,E2,E3,E1dir,E2dir,E3dir] = EigCalc3d(Exx,Eyy,Ezz,Exy,Exz,Eyz)
%
% Arguments: (input)
% Sxx,Syy,Szz
% Sxy,Sxz,Syz       - The stress tensor components at a point, each can be a
%                    column vector.
%
% Arguments: (output)
% S1,S2,S3       	- Principal stress component magnitudes Sigma 1,2 and
%                    Sigma 3.
%
% S1dir,S2dir,S3dir - Principal stress directions (direction cosines). Each
%                    will be a n*3 column vector [CosAx,CosAy,CosAz] of this
%                    direction. 
%
% Example usage 1:
%
% %Calculating directions for a 3D stress tensor
% Sxx=0.2; Syy=-1.5; Szz=0;
% Sxy=1;   Sxz=0;    Syz=0;
% [S1,S2,S3,S1dir,S2dir,S3dir] = EigCalc3d(Sxx,Syy,Szz,Sxy,Sxz,Syz)
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Preallocating array
n = numel(Sxx);
tensor=zeros (n*3,3);

%Filling this array with 3x3 tensors, every 3 rows is each point
%Indexing way to accumulate mats
N=3;
tensor(1:N:end,1) = Sxx(1:1:end,:);
tensor(1:N:end,2) = Sxy(1:1:end,:);
tensor(1:N:end,3) = Sxz(1:1:end,:);
tensor(2:N:end,1) = Sxy(1:1:end,:);
tensor(2:N:end,2) = Syy(1:1:end,:);
tensor(2:N:end,3) = Syz(1:1:end,:);
tensor(3:N:end,1) = Sxz(1:1:end,:);
tensor(3:N:end,2) = Syz(1:1:end,:);
tensor(3:N:end,3) = Szz(1:1:end,:);

%Eig can't handle nan's so we turn these to 0's and put the calculated s1s2s3 to nans after 
NanFlag = isnan(tensor);
tensor(NanFlag)=0;

%Calculating the eigen vectors and principal values (diagonal)
%Preallocating array
V=zeros (n*3,3); %Eigen vectors - http://eqseis.geosc.psu.edu/~cammon/HTML/UsingMATLAB/PDF/ML2%20Eigenvalues.pdf
D=zeros (n*3,3); %Eigen values -  also see Pollard 2005 addition resources MATLAB introduction WordDoc

for J=(1:3:n*3) 
    [V(J:J+2,:),D(J:J+2,:)] = eig(tensor(J:J+2,:));
    V(J:J+2,:)=V(J:J+2,:)'; %Flipping the direction cosines as its easier to extract these like this. 
end


%Putting anywhere where there were nans in the tensor to nan
D(NanFlag)=nan;
V(NanFlag)=nan;

%Now getting a col vec of S3 S2 S1 which is for each point
%J=(1:3:size(D(:,1)));
J = 1:3:length(D(:,1));
for i = 1:n 
    [B(J(:,i):J(:,i)+2,:),~]=sort(diag(D(J(:,i):J(:,i)+2,:)));
end  

%B output is S1S2S3 in a single column list, s1(a),s2(a),s3(a),s1(b),s2(b) etc. S1 S2 and S3 are then extracted to thier own arrays ready for 
%export from the function.
N=3;
S3 = B(1:N:end,:);
S2 = B(2:N:end,:);
S1 = B(3:N:end,:);

%Also getting the direction cosines in the same manner. (note that X is first col, Y is second etc) 
S3dir = V(1:N:end,:);
S2dir = V(2:N:end,:);
S1dir = V(3:N:end,:);

% How to draw fig with scaled vectors (have to do outside this func where you have a XYZ)
% figure;a = quiver3(X(:),Y(:),Z(:),S1dirt(:,1),S1dirt(:,2),S1dirt(:,3));
% hold on
% a2 = quiver3(X(:),Y(:),Z(:),-S1dirt(:,1),-S1dirt(:,2),-S1dirt(:,3)); %other dir
% b = quiver3(X(:),Y(:),Z(:),S2dirt(:,1),S2dirt(:,2),S2dirt(:,3));
% b2 = quiver3(X(:),Y(:),Z(:),-S2dirt(:,1),-S2dirt(:,2),-S2dirt(:,3)); %other dir
% c = quiver3(X(:),Y(:),Z(:),S3dirt(:,1),S3dirt(:,2),S3dirt(:,3));
% c2 = quiver3(X(:),Y(:),Z(:),-S3dirt(:,1),-S3dirt(:,2),-S3dirt(:,3)); %other dir
% title({'\fontsize{14}Principal directions','\fontsize{8}Red S1, Blue S2, Green S3, scaled by magnitude'})
% axis equal ;ax.projection = 'perspective';
% Lw=2;%Line width prop, default is .5
% a.Color = 'red'  ;a.ShowArrowHead='off';a.LineWidth = Lw; a2.Color = 'red'  ;a2.ShowArrowHead='off';a2.LineWidth = Lw;
% b.Color = 'green';b.ShowArrowHead='off';b.LineWidth = Lw; b2.Color = 'green';b2.ShowArrowHead='off';b2.LineWidth = Lw;
% c.Color = 'blue' ;c.ShowArrowHead='off';c.LineWidth = Lw; c2.Color = 'blue' ;c2.ShowArrowHead='off';c2.LineWidth = Lw;
% %set(a,'Visible','off');set(a2,'Visible','off') %To turn off certain things do this

end

