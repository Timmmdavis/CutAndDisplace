function [S1,S2,S1dir,S2dir]=EigCalc2d(Sxx,Syy,Sxy)
% Calculates the 2-D principal stress/strain magnitudes and directions from input tensors
% Input parameters
% Sxx = sigma xx
% Syy = sigma yy
% Sxy = sigma xy
% Output parameter (tension positive convention for inputs)
% S1 = most tensile principal stress
% S2 = most comp principal stress
% S1dir = XY components of the direction of the principal stress S1
% S2dir = see above
% Example 
%Calclating 2d stress eigenValues
%[S1,S2,S1dir,S2dir]=EigCalc2d(Sxx,Syy,Sxy);
%Calclating 2d strain eigenValues
%[E1,E2,E1dir,E2dir]=EigCalc2d(Exx,Eyy,Exy);

%   Copyright 2017, Tim Davis, The University of Aberdeen

%Preallocating array
n = numel(Sxx);
tensor=zeros (n*2,2);
%Filling this array with 3x3 tensors, every 3 rows is each point
J=(1:2:n*2);
for i = 1:n 
    
    tensor(J(:,i):J(:,i)+1,:)=[Sxx(i), Sxy(i);
                               Sxy(i), Syy(i)]; 
                 
end

%Eig can't handle nan's so we turn these to 0's and put the calculated s1s2 to nans after 
NanFlag = isnan(tensor);
tensor(NanFlag)=0;
NanFlag = isinf(tensor);
tensor(NanFlag)=0;

%Calculating the eigen vectors and principal values (diagonal)
%Preallocating array
V=zeros (n*2,2); %Eigen vectors - http://eqseis.geosc.psu.edu/~cammon/HTML/UsingMATLAB/PDF/ML2%20Eigenvalues.pdf
D=zeros (n*2,2); %Eigen values -  also see Pollard 2005 addition resources MATLAB introduction WordDoc
for J=(1:2:n*2); 
    [V(J:J+1,:),D(J:J+1,:)] = eig(tensor(J:J+1,:));
    V(J:J+1,:)=V(J:J+1,:)'; %Flipping the direction cosines as its easier to extract these like this. 
end

%Putting anywhere where there were nans in the tensor to nan
D(NanFlag)=nan;
V(NanFlag)=nan;

%Now getting a col vec of S3 S2 S1 which is for each point
J=(1:2:size(D(:,1)));
for i = 1:n 
    [B(J(:,i):J(:,i)+1,:),I(J(:,i):J(:,i)+1,:)]=sort(diag(D(J(:,i):J(:,i)+1,:)));
end  

%B output is S1S2S3 in a single column list, s1(a),s2(a),s3(a),s1(b),s2(b) etc. S1 S2 and S3 are then extracted to thier own arrays ready for 
%export from the function.
N=2;
S2 = B(1:N:end,:);
S1 = B(2:N:end,:);

%Also getting the direction cosines in the same manner. (note that X is first col, Y is second etc) 
S2dir = V(1:N:end,:);
S1dir = V(2:N:end,:);

%The other way of doing this would be:
% S1 = (Sxx+Syy)/2 + sqrt( ((Sxx-Syy)/2).^2 + Sxy.^2);
% y1=(1/2)*atan2(2*Sxy, Sxx-Syy); % pg2 222 pollard eq6.73
% S2 = (Sxx+Syy)/2 - sqrt( ((Sxx-Syy)/2).^2 + Sxy.^2); (just add 1/4 pi r
% to y1 to get y2). 

