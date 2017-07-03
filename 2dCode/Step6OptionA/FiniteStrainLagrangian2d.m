function [exx,eyy,exy,ExxInfErrorPerc,EyyInfErrorPerc,ExyInfErrorPerc]=FiniteStrainLagrangian2d(Ux,Uy,X,Y)
% Calculate the Lagrangian strain from input grid with displacement vectors (displacement from original position). 
%
% (3D)  E = STRAIN(Ux,Uy,X,Y)
%
% inputs,
%   Ux,Uy: The displacement vector in the
%           x and y direction
%   X,Y:   The grid of points, can be col vecs of sparsely distributed
%          points or a nxn square/rectangular array
%
% outputs,
%   tensors: the 2-D Lagrangian strain tensors, [exx,eyy,exy]; %FINITE
%            if you want infinitesimal ones these are computed also and
%            could be exported
%   error:   the % error an infinitesimal strain calculation will introduce
%            values above 100% have been set to 100.
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

[ Uniform ] = UniformGridCheck2d( X,Y );
if Uniform==1;

        %%Calculating spacing
    spacingx=X(1,1)-X(1,2);
    spacingy=Y(1,1)-Y(2,1);
    
        %%If Ux and UY were put in as cols we fix this here
    if iscolumn(Ux)
        Ux=reshape(Ux,size(X));
    end
    if iscolumn(Uy)
        Uy=reshape(Uy,size(X));
    end

% %     % Calculating displacement gradients, This is not deformation gradient
% %     % F, Martel lecture 14 GG303 page 6. See 15 page 3+4
    [Uxx,Uxy] = gradient(Ux,abs(spacingx)); %changed from D.Kroons to match Steve Martels lecture notes
    [Uyx,Uyy] = gradient(Uy,abs(spacingy));
    
else %could use interpolation instead of this method but that can give dodgy edge gradients. Neither are perfect
    
    X=X(:);
    Y=Y(:);
    Ux=Ux(:);
    Uy=Uy(:);
    Filler=zeros(numel(X),6);

    disp('calculating gradient matrices for non uniform grid')
    for i=1:numel(X)  

    %only getting points closest to first 
    A = [X(i),Y(i)]; 
    B = [X(:),Y(:)];
    dist=bsxfun(@hypot,B(:,1)-A(1),B(:,2)-A(2)); %total distance from point A
    B = [X(:),Y(:),Ux(:),Uy(:)];    %arranging vec
    B = sortrows([dist,B]);         %sorting on distance, smallest first
    %B is now a vector of [distance from A, X vals, Y vals,Ux,Uy]
    %First row can be ignored as it is the same point as A
    closestPoints=B(2:21,2:5); %the 20 closest points to the one we are calculating strain on

    %Extracting vars
    C=closestPoints(:,1:2); %[X,Y] coordinates
    d=closestPoints(:,3:4); %[Ux,Uy]


    % get velocity data and reformat into column vector d
    d = d'; % velocities by row then column
    d = d(:); % make column vector
    % Every 1st row is Ux for the first point and the second is Uy 


    % set up coefficient matrix (Nestor Style!)
    G2 = [1, 0, C(1,2),  C(1,1),    0,      0     ;...
         0, 1,  0,       0,         C(1,2), C(1,1);...
         1, 0,  C(2,2), C(2,1),     0,      0     ;...
         0, 1,  0,      0,          C(2,2), C(2,1);...
         1, 0,  C(3,2), C(3,1),     0,      0     ;...
         0, 1,  0,      0,          C(3,2), C(3,1);
         1, 0,  C(4,2), C(4,1),     0,      0     ;...
         0, 1,  0,      0,          C(4,2), C(4,1); 
         1, 0,  C(5,2), C(5,1),    0,      0     ;...
         0, 1,  0,       0,         C(5,2), C(5,1);...
         1, 0,  C(6,2), C(6,1),     0,      0     ;...
         0, 1,  0,      0,          C(6,2), C(6,1);...
         1, 0,  C(7,2), C(7,1),     0,      0     ;...
         0, 1,  0,      0,          C(7,2), C(7,1);
         1, 0,  C(8,2), C(8,1),     0,      0     ;...
         0, 1,  0,      0,          C(8,2), C(8,1); 
         1, 0,  C(9,2),  C(9,1),    0,      0     ;...
         0, 1,  0,       0,         C(9,2), C(9,1);...
         1, 0,  C(10,2), C(10,1),     0,      0     ;...
         0, 1,  0,      0,          C(10,2), C(10,1);...
         1, 0,  C(11,2), C(11,1),     0,      0     ;...
         0, 1,  0,      0,          C(11,2), C(11,1);
         1, 0,  C(12,2), C(12,1),     0,      0     ;...
         0, 1,  0,      0,          C(12,2), C(12,1); 
         1, 0,  C(13,2), C(13,1),    0,      0     ;...
         0, 1,  0,       0,         C(13,2), C(13,1);...
         1, 0,  C(14,2), C(14,1),     0,      0     ;...
         0, 1,  0,      0,          C(14,2), C(14,1);...
         1, 0,  C(15,2), C(15,1),     0,      0     ;...
         0, 1,  0,      0,          C(15,2), C(15,1);
         1, 0,  C(16,2), C(16,1),     0,      0     ;...
         0, 1,  0,      0,          C(16,2), C(16,1); 
         1, 0,  C(17,2),  C(17,1),    0,      0     ;...
         0, 1,  0,       0,         C(17,2), C(17,1);...
         1, 0,  C(18,2), C(18,1),     0,      0     ;...
         0, 1,  0,      0,          C(18,2), C(18,1);...
         1, 0,  C(19,2), C(19,1),     0,      0     ;...
         0, 1,  0,      0,          C(19,2), C(19,1);
         1, 0,  C(20,2), C(20,1),     0,      0     ;...
         0, 1,  0,      0,          C(20,2), C(20,1)];

    % solve for model paramaters (tx, ty, Uxy, Uxx, Uyy, Uyx)
    m=G2\d;

    Filler(i,:)=m';  %filling each row with the calculated unknowns

    %simple progress
    (i/numel(X))*100 %very simple progress prints to cmd and counts to 100

    end

    Uxy=Filler(:,3);
    Uxx=Filler(:,4);
    Uyy=Filler(:,5);
    Uyx=Filler(:,6);

end %end the calculation of the gradient matrices

    %Now reshaping these displacement gradients into column vectors (square
    %if we used the grid option).
    [Uxx] = Uxx(:);
    [Uxy] = Uxy(:);
    [Uyx] = Uyx(:);
    [Uyy] = Uyy(:);

    %Creating a list of displacement tensors, each 3x3 row is a different point. 
    n=numel(Ux);
    J=(1:2:n*2);
    for i = 1:n 
    Ugrad(J(:,i):J(:,i)+1,:)=[Uxx(i), Uxy(i);
                              Uyx(i), Uyy(i)];           
    end
    
    %Now calculating the transformation matrix F, note this removes inf
    %values that can appear at the model edges that have strange gradients. 
    %The eig function returns an error when it meets inf values.  
    %http://www.files.ethz.ch/structuralgeology/sms/NumModRocDef/Strain_Tensors.pdf
    %eq 3
    for i = 1:n 
    F(J(:,i):J(:,i)+1,:)=Ugrad(J(:,i):J(:,i)+1,:)+[1 0 ;0 1]; %F=inv(Finv);
    Ft(J(:,i):J(:,i)+1,:)=F(J(:,i):J(:,i)+1,:)';
    end
    
    bad = isinf(F);
    F(bad) = 0;
    bad = isinf(Ft);
    Ft(bad) = 0;
    
    %Calculating Lagrangian strain. Using the formula from Martel
    %'GG303/Lec15_2015 page 6-9
    E=zeros(n*2,2);
    for i = 1:n 
    E(J(:,i):J(:,i)+1,:)=(1/2)*((Ft(J(:,i):J(:,i)+1,:)*F(J(:,i):J(:,i)+1,:))-[1 0;0 1]);
    %Turn on to calculating infinitesimal strain. Using the formula from Martel
    %'GG303/Lec15_2015 page 9
    e(J(:,i):J(:,i)+1,:)=(1/2)*(Ugrad(J(:,i):J(:,i)+1,:)+(Ugrad(J(:,i):J(:,i)+1,:)'));
    end
   
    %extracting the individual tensors as column vectors 
    N=2;
    exx= E(1:N:end,1);
    exy= E(1:N:end,2);
    eyx= E(2:N:end,1); %same as exy
    eyy= E(2:N:end,2);

    
    %Using calculated infinitesimal strain to calcluate error associated
    %with infinitesimal calcluation compared to finite. 
    %See D Pollard  2005 Chapter 5.5.3
    exxinf= e(1:N:end,1);
    exyinf= e(1:N:end,2);
    eyyinf= e(2:N:end,2);
    ExxInfErrorPerc=abs(((exx-exxinf)./exx)*100);
    EyyInfErrorPerc=abs(((eyy-eyyinf)./eyy)*100);
    ExyInfErrorPerc=abs(((exy-exyinf)./exy)*100);
    
    %Making sure max is 100% anymore and we put to Nan. Can get crazy
    %values across the fault where we don’t really care
    bad=abs(ExxInfErrorPerc)>100;
    ExxInfErrorPerc(bad)=100;
    bad=abs(EyyInfErrorPerc)>100;
    EyyInfErrorPerc(bad)=100;
    bad=abs(ExyInfErrorPerc)>100;
    ExyInfErrorPerc(bad)=100;

    
    %Turn off options below if you don’t want figures
    
    %Drawing error maps
    X=X(:);
    Y=Y(:);
    figure;subplot(1,3,1),scatter(X,Y,75,ExxInfErrorPerc,'filled', 's');h = colorbar;title(h,'100+')
    xlabel('x'); ylabel('z'); axis('equal'); title('ExxInfinitesimalError%'); axis([min(X),max(X),min(Y),max(Y)])
    subplot(1,3,2),scatter(X,Y,75,EyyInfErrorPerc,'filled', 's');h = colorbar;title(h,'100+')
    xlabel('x'); ylabel('z'); axis('equal'); title('EyyInfinitesimalError%'); axis([min(X),max(X),min(Y),max(Y)])
    subplot(1,3,3),scatter(X,Y,75,ExyInfErrorPerc,'filled', 's');h = colorbar;title(h,'100+')
    xlabel('x'); ylabel('z'); axis('equal'); title('ExyInfinitesimalError%'); axis([min(X),max(X),min(Y),max(Y)])

    %Drawing in Finite strain
    ClipAbove=0.25; %remove values above this limit, good for points across fault surfaces
    bad=abs(exx)>ClipAbove;
    exx(bad)=nan;
    bad=abs(eyy)>ClipAbove;
    eyy(bad)=nan;
    bad=abs(exy)>ClipAbove;
    exy(bad)=nan;
    figure;subplot(1,3,1),scatter(X,Y,75,exx,'filled', 's');colorbar;
    xlabel('x'); ylabel('z'); axis('equal'); title('Exx'); axis([min(X),max(X),min(Y),max(Y)])
    subplot(1,3,2),scatter(X,Y,75,eyy,'filled', 's');colorbar;
    xlabel('x'); ylabel('z'); axis('equal'); title('Eyy'); axis([min(X),max(X),min(Y),max(Y)])
    subplot(1,3,3),scatter(X,Y,75,exy,'filled', 's');colorbar;
    xlabel('x'); ylabel('z'); axis('equal'); title('Exy'); axis([min(X),max(X),min(Y),max(Y)])



 
