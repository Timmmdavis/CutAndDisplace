function [S1,S2,S1dir,S2dir]=EigCalc2dSpeed(Sxx,Syy,Sxy)
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
tensor=zeros (2,2,n);

%Filling this array with 2x2 tensors, each 3rd dim is each point
%Accumulated using indexing
%[Sxx,Sxy]
%[Syx,Syy]
tensor(1,1,1:1:end) = Sxx(1:1:end,:);
tensor(1,2,1:1:end) = Sxy(1:1:end,:);
tensor(2,1,1:1:end) = Sxy(1:1:end,:);
tensor(2,2,1:1:end) = Syy(1:1:end,:);

%Eig can't handle nan's so we turn these to 0's and put the calculated s1s2 to nans after 
NanFlag = isnan(tensor);
tensor(NanFlag)=0;
NanFlag = isinf(tensor);
tensor(NanFlag)=0;

%Do the calculation
D = eig2(tensor); %see base of script and functions it calls
%Sort results
D = sortrows(D);
%Transpose
D=D';
%Extract
S2 = D(:,1);
S1 = D(:,2);

%Create vars for loop (preallocate)
S1dir=zeros(2,n);
%Create Identity matrix
Ident = repmat(eye(2),1,1,n);

for i=1:n
    %Equation we want to solve
    %[A-(Lambda * Identity)]v=0
    
    %Grab the current bit in the loop
    I=Ident(:,:,i);
    A=tensor(:,:,i);
    
    %[Lambda * Identity] in eigenvector equation
    One=S1(i)*I;
    
    %[A-(Lambda * Identity)]
    vect1=A-One;
    
    %[A-(Lambda * Identity)]v=0. Finding vector [v]
    %https://de.mathworks.com/help/matlab/ref/decomposition.html
    [Q1,~] = qr(vect1,0);

    %Extracting results
    S1dir(:,i)=(Q1(:,2));
    
end

S1dir=S1dir';

%We know S2 is perpendicular to S1. 
S2dir=[S1dir(2,:),-S1dir(1,:)];

%The other way of doing this would be:
% S1 = (Sxx+Syy)/2 + sqrt( ((Sxx-Syy)/2).^2 + Sxy.^2);
% y1=(1/2)*atan2(2*Sxy, Sxx-Syy); % pg2 222 pollard eq6.73
% S2 = (Sxx+Syy)/2 - sqrt( ((Sxx-Syy)/2).^2 + Sxy.^2); (just add 1/4 pi r
% to y1 to get y2). 

