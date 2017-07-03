function [exx,eyy,ezz,exy,exz,eyz,ExxInfErrorPerc,EyyInfErrorPerc,EzzInfErrorPerc,ExyInfErrorPerc,ExzInfErrorPerc,EyzInfErrorPerc]=FiniteStrainLagrangian3d(Ux,Uy,Uz,X,Y,Z)
% Calculate the Lagrangian strain from input grid with displacement vectors (displacement from original position). 
%
% (3D)  E = STRAIN(Ux,Uy,X,Y)
%
% inputs,
%   Ux,Uy,Uz: The displacement vector in the
%           x , y and z direction
%   X,Y,Z:   The grid of points, can be col vecs of sparsely distributed
%          points or a nxn square/rectangular array
%
% outputs,
%   tensors: the 3-D Lagrangian strain tensors, [exx,eyy,ezz,exy,exz,eyz]; %FINITE
%            if you want infinitesimal ones these are computed also and
%            could be exported
%   error:   the % error an infinitesimal strain calculation will introduce
%
% Sources used:
%   For the gridded function (much faster) the idea to use gradient is from
%   a function written by D.Kroon University of Twente (February 2009)
%   found on the MATLAB file exchange
%   For the non uniform data the method used is that described in 'Cardozo,
%   N. and Allmendinger, R.W., 2009. SSPX: A program to compute strain from
%   displacement/velocity data. Computers & Geosciences, 35(6),
%   pp.1343-1357.I would not have been able to work this out so quickly
%   though without the 'calcstrain.m file written by Phil Resor for
%   infinitesimal strain on January 7, 2013
%   The lecture notes of Steve Martel at SOESST where particularly clear
%   and helpful for the actual strain calculations, as were Marcel
%   Frehner's notes on strain (University of Zurich)
%   Ryan Shackleton helped in my understanding of this topic.

%   Copyright 2017, Tim Davis, The University of Aberdeen
%These loops check if the data is gridded, if so it will use the gridded
%gradient matrix calculation. Flagx and Flagy are created. If any of these are 1
%the data is not a uniform grid. The precision is set to 8 (I round) so
%really small grids will suffer here.
%Grid data check X
[ Uniform ] = UniformGridCheck3d( X,Y,Z );

if Uniform==1

        %%Calculating spacing
    spacingx=X(1,1,1)-X(1,2,1);
    spacingy=Y(1,1,1)-Y(2,1,1);
    spacingz=Z(1,1,1)-Z(1,1,2);
    

        %%If Ux and UY were put in as cols we fix this here
    if iscolumn(Ux)
        Ux=reshape(Ux,size(X));
    end
    if iscolumn(Uy)
        Uy=reshape(Uy,size(X));
    end
    if iscolumn(Uz)
        Uz=reshape(Uz,size(X));
    end

% %     % Calculating displacement gradients, This is not deformation gradient
% %     % F, Martel lecture 14 GG303 page 6. See 15 page 3+4
    [Uxx,Uxy,Uxz] = gradient(Ux,spacingx); %changed from D.Kroons to match Steve Martels lecture notes
    [Uyx,Uyy,Uyz] = gradient(Uy,spacingy);
    [Uzx,Uzy,Uzz] = gradient(Uz,spacingz);
    

    
else %could use interpolation instead of this method but that can give dodgy edge gradients. Neither are perfect
    
    X=X(:);
    Y=Y(:);
    Z=Z(:);
    Ux=Ux(:);
    Uy=Uy(:);
    Uz=Uz(:);
    Filler=zeros(numel(X),12);

    disp('calculating gradient matrices for non uniform grid')
    for i=1:numel(X)  

    %only getting points closest to first 
    A = [X(i),Y(i),Z(i)]; 
    B = [X(:),Y(:),Z(:)];
    x2=(B(:,1)-A(1)).^2; %length x
    y2=(B(:,2)-A(2)).^2; %length y
    z2=(B(:,3)-A(3)).^2;
    dist=sqrt(abs(x2)+abs(y2)+abs(z2)); %total distance from point A
    B = [X(:),Y(:),Z(:),Ux(:),Uy(:),Uz(:)];    %arranging vec
    B = sortrows([dist,B]);         %sorting on distance, smallest first
    %B is now a vector of [distance from A, X vals, Y vals, Zvals, Ux,Uy,Uz]
    %First row can be ignored as it is the same point as A
    closestPoints=B(2:21,2:7); %the 20 closest points to the one we are calculating strain on

    %Extracting vars
    C=closestPoints(:,1:3); %[X,Y,Z] coordinates
    d=closestPoints(:,4:6); %[Ux,Uy,Uz]


    % get velocity data and reformat into column vector d
    d = d'; % velocities by row then column
    d = d(:); % make column vector
    % Every 1st row is Ux for the first point and the second is Uy third is Uz 


    % set up coefficient matrix (Nestor Style!)
    G2 = [1, 0, 0,  C(1,3),  C(1,2),  C(1,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(1,3), C(1,2), C(1,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(1,3), C(1,2), C(1,1);...
          1, 0, 0,  C(2,3),  C(2,2),  C(2,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(2,3), C(2,2), C(2,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(2,3), C(2,2), C(2,1);...
          1, 0, 0,  C(3,3),  C(3,2),  C(3,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(3,3), C(3,2), C(3,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(3,3), C(3,2), C(3,1);...
          1, 0, 0,  C(4,3),  C(4,2),  C(4,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(4,3), C(4,2), C(4,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(4,3), C(4,2), C(4,1);...
          1, 0, 0,  C(5,3),  C(5,2),  C(5,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(5,3), C(5,2), C(5,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(5,3), C(5,2), C(5,1);...
          1, 0, 0,  C(6,3),  C(6,2),  C(6,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(6,3), C(6,2), C(6,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(6,3), C(6,2), C(6,1);...
          1, 0, 0,  C(7,3),  C(7,2),  C(7,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(7,3), C(7,2), C(7,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(7,3), C(7,2), C(7,1);... 
          1, 0, 0,  C(8,3),  C(8,2),  C(8,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(8,3), C(8,2), C(8,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(8,3), C(8,2), C(8,1);...
          1, 0, 0,  C(9,3),  C(9,2),  C(9,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(9,3), C(9,2), C(9,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(9,3), C(9,2), C(9,1);...       
          1, 0, 0,  C(10,3),  C(10,2),  C(10,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(10,3), C(10,2), C(10,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(10,3), C(10,2), C(10,1);...             
          1, 0, 0,  C(11,3),  C(11,2),  C(11,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(11,3), C(11,2), C(11,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(11,3), C(11,2), C(11,1);...
          1, 0, 0,  C(12,3),  C(12,2),  C(12,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(12,3), C(12,2), C(12,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(12,3), C(12,2), C(12,1);...
          1, 0, 0,  C(13,3),  C(13,2),  C(13,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(13,3), C(13,2), C(13,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(13,3), C(13,2), C(13,1);...
          1, 0, 0,  C(14,3),  C(14,2),  C(14,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(14,3), C(14,2), C(14,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(14,3), C(14,2), C(14,1);...
          1, 0, 0,  C(15,3),  C(15,2),  C(15,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(15,3), C(15,2), C(15,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(15,3), C(15,2), C(15,1);...
          1, 0, 0,  C(16,3),  C(16,2),  C(16,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(16,3), C(16,2), C(16,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(16,3), C(16,2), C(16,1);...
          1, 0, 0,  C(17,3),  C(17,2),  C(17,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(17,3), C(17,2), C(17,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(17,3), C(17,2), C(17,1);... 
          1, 0, 0,  C(18,3),  C(18,2),  C(18,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(18,3), C(18,2), C(18,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(18,3), C(18,2), C(18,1);...
          1, 0, 0,  C(19,3),  C(19,2),  C(19,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(19,3), C(19,2), C(19,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(19,3), C(19,2), C(19,1);...          
          1, 0, 0,  C(20,3),  C(20,2),  C(20,1),    0,      0,      0,       0,      0,      0,    ;...
          0, 1, 0,  0,       0,       0,         C(20,3), C(20,2), C(20,1),  0,      0,      0,    ;...
          0, 0, 1,  0,       0,       0,         0,      0,      0,       C(20,3), C(20,2), C(20,1)];
      
    % solve for model paramaters (tx, ty, Uxy, Uxx, Uyy, Uyx)
    m=G2\d;

    Filler(i,:)=m';  %filling each row with the calculated unknowns

    %simple progress
    (i/numel(X))*100 %very simple progress prints to cmd and counts to 100

    end

    %extracting, matrix is back to front compared to normal notation but
    %this works fine
    Uxx=-Filler(:,6);
    Uyy=-Filler(:,8);   
    Uzz=-Filler(:,10);
    Uxy=-Filler(:,5); 
    Uxz=-Filler(:,4);    
    Uyz=-Filler(:,7);
    Uyx=-Filler(:,9);
    Uzx=-Filler(:,12);
    Uzy=-Filler(:,11);


end %end the calculation of the gradient matrices

    %Now reshaping these displacement gradients into column vectors. 
    Uxx=Uxx(:);
    Uxy=Uxy(:);
    Uxz=Uxz(:);
    Uyx=Uyx(:);
    Uyy=Uyy(:);
    Uyz=Uyz(:);
    Uzx=Uzx(:);
    Uzy=Uzy(:);
    Uzz=Uzz(:);
    
    n=numel(Ux);
   %Creating a list of displacement tensors, each 3x3 row is a different point. 
    J=(1:3:n*3);
    for i = 1:n 
    Ugrad(J(:,i):J(:,i)+2,:)=[Uxx(i), Uxy(i), Uxz(i);
                               Uyx(i), Uyy(i), Uyz(i);
                               Uzx(i), Uzy(i), Uzz(i)];           
    end
    
    
    %Now calculating the transformation matrix F, note this removes inf
    %values that can appear at the model edges that have strange gradients. 
    %The eig function returns an error when it meets inf values.  
    %http://www.files.ethz.ch/structuralgeology/sms/NumModRocDef/Strain_Tensors.pdf
    %eq 3
    for i = 1:n 
    F(J(:,i):J(:,i)+2,:)=Ugrad(J(:,i):J(:,i)+2,:)+[1 0 0;0 1 0;0 0 1]; %F=inv(Finv);
    Ft(J(:,i):J(:,i)+2,:)=F(J(:,i):J(:,i)+2,:)';
    end
    
    bad = isinf(F);
    F(bad) = 0;
    bad = isinf(Ft);
    Ft(bad) = 0;
    
    %Calculating Lagrangian strain. Using the formula from Martel
    %'GG303/Lec15_2015 page 6-9
    E=zeros(n*3,3);
    for i = 1:n 
    E(J(:,i):J(:,i)+2,:)=(1/2)*((Ft(J(:,i):J(:,i)+2,:)*F(J(:,i):J(:,i)+2,:))-[1 0 0;0 1 0;0 0 1]);
    %Turn on to calculating infinitesimal strain. Using the formula from Martel
    %'GG303/Lec15_2015 page 9
    e(J(:,i):J(:,i)+2,:)=(1/2)*(Ugrad(J(:,i):J(:,i)+2,:)+(Ugrad(J(:,i):J(:,i)+2,:)'));
    end
   
    
    %extracting the individual tensors as column vectors 
    N=3;
    exx= E(1:N:end,1);
    exy= E(1:N:end,2);
    exz= E(1:N:end,3);
    eyx= E(2:N:end,1); %same as exy
    eyy= E(2:N:end,2);
    eyz= E(2:N:end,3);
    ezx= E(3:N:end,1); %same as ezx
    ezy= E(3:N:end,2); %same as ezy
    ezz= E(3:N:end,3);
    
    %Using calculated infinitesimal strain to calcluate error associated
    %with infinitesimal calcluation compared to finite. 
    %See D Pollard  2005 Chapter 5.5.3
    exxinf= e(1:N:end,1);
    exyinf= e(1:N:end,2);
    exzinf= e(1:N:end,3);
    eyyinf= e(2:N:end,2);
    eyzinf= e(2:N:end,3);
    ezzinf= e(3:N:end,3);
    ExxInfErrorPerc=((exx-exxinf)./exx)*100;
    EyyInfErrorPerc=((eyy-eyyinf)./eyy)*100;
    EzzInfErrorPerc=((ezz-ezzinf)./ezz)*100;
    ExyInfErrorPerc=((exy-exyinf)./exy)*100;
    ExzInfErrorPerc=((exz-exzinf)./exz)*100;
    EyzInfErrorPerc=((eyz-eyzinf)./eyz)*100;
    
    
    
    




 
