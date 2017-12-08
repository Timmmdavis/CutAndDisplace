function [exx,eyy,ezz,exy,exz,eyz,ExxIErr,EyyIErr,EzzIErr,ExyIErr,ExzIErr,EyxIErr]=FiniteStrainLagrangian3d(Ux,Uy,Uz,X,Y,Z)
% FiniteStrainLagrangian3d: Calculate the Lagrangian strain from input grid 
%                   with displacement vectors (displacement from original
%                   position). Note the data does not need to be gridded.
%                   Non gridded data calculates strain based on the 20
%                   closest points. This also supplies the error in %
%                   introduced due to the assumption of infinitesimal
%                   strain.
%
%                   Sources used:
%                   For the gridded function (much faster) the idea to use
%                   gradient is from a function written by D.Kroon
%                   University of Twente (February 2009) on the MATLAB file
%                   exchange.
%                   For the non uniform data the method used is that
%                   described in 'Cardozo, N. and Allmendinger, R.W., 2009.
%                   SSPX: A program to compute strain from
%                   displacement/velocity data. Computers & Geosciences,
%                   35(6), pp.1343-1357.
%                   I would not have been able to work this out so quickly
%                   though without the 'calcstrain.m' file written by Phil
%                   Resor for infinitesimal strain on January 7, 2013.
%                   The lecture notes of Steve Martel at SOESST where
%                   particularly clear and helpful for the actual strain
%                   calculations, as were Marcel Frehner's notes on strain
%                   (University of Zurich).
%                   Error calculations from Pollard and Fletcher, 2005.
%                   Ryan Shackleton helped in my understanding of this
%                   topic.
%   
% usage #1: 
% [exx,eyy,ezz,exy,exz,eyz,ExxIErr,EyyIErr,EzzIErr,ExyIErr,ExzIErr,EyxIErr]=FiniteStrainLagrangian3d(Ux,Uy,Uz,X,Y,Z)
%
% usage #2: If the infinitesimal strain errors dont concern you. 
% [exx,eyy,ezz,exy,exz,eyz]=FiniteStrainLagrangian3d(Ux,Uy,Uz,X,Y,Z)
%
% Arguments: (input)
% Ux,Uy,Uz          - Displacement components in X Y and Z at points defined
%                    in the vectors X Y and Z.
%
% X,Y,Z             - The location of the points in 3D space.
%
% Arguments: (output)
% exx,eyy,ezz
% exy,exz,eyz       - The 3-D Lagrangian finite strain tensors,
%                    if you want infinitesimal ones these are computed also
%                    and could be exported.
%
% ExxIErr,EyyIErr
% EzzIErr,ExyIErr
% ExzIErr,EyzIErr   - Error in percent if we assume infinitesimal strain. 
%
%
% Example usage 1: Gridded data:
% 
%  [X,Y,Z]=meshgrid(-2:0.3:2, -2:0.3:2, -2:0.3:2);
%  dimx = length(X(:,:,1)); 
%  dimy = length(X(:,1,:));
%  dimz = length(X(1,:,:));
%  [~,~,r] = cart2sph(X(:),Y(:),Z(:));
%  Ux=(r(:).^2).*X(:);
%  Uy=(r(:).^2).*Y(:);
%  Uz=(r(:).^2).*Z(:);
%  %Only grabbing the normal components
%  [exx,eyy,ezz]=FiniteStrainLagrangian3d(Ux,Uy,Uz,X,Y,Z);
%  [exx,eyy,ezz]=ReshapeData3d( dimx,dimy,dimz,exx,eyy,ezz);
%  IsoContoursPrincipalStressPercentiles( exx,eyy,ezz,X,Y,Z);
%  %Removing titles as they are not correct for this case:
%  delete(findall(findall(gcf,'Type','axe'),'Type','text'))
% 
% Example usage 2: Non gridded data:
%
% Density=5000;
% n=4;
% xmv=-2; ymv=-2;zmv=-2;
% X=(rand(1,Density)*n)'+xmv;
% Y=(rand(1,Density)*n)'+ymv;
% Z=(rand(1,Density)*n)'+zmv;
% [~,~,r] = cart2sph(X(:),Y(:),Z(:));
% Ux=(r.^2).*X;
% Uy=(r.^2).*Y;
% Uz=(r.^2).*Z;
% %Only grabbing the normal components
% [exx,eyy,ezz]=FiniteStrainLagrangian3d(Ux,Uy,Uz,X,Y,Z);
% %Now interpolating this so we can draw it:
% [Xg,Yg,Zg]=meshgrid(-2:0.3:2, -2:0.3:2, -2:0.3:2);
% dimx = length(X(:,:,1)); 
% dimy = length(X(:,1,:));
% dimz = length(X(1,:,:));
% Fexx = scatteredInterpolant(X,Y,Z,exx);
% Vexxg = Fexx(Xg,Yg,Zg);
% Feyy = scatteredInterpolant(X,Y,Z,eyy);
% Veyyg = Feyy(Xg,Yg,Zg);
% Fezz = scatteredInterpolant(X,Y,Z,ezz);
% Vezzg = Fezz(Xg,Yg,Zg);
% IsoContoursPrincipalStressPercentiles( Vexxg,Veyyg,Vezzg,Xg,Yg,Zg);
% %Removing titles as they are not correct for this case:
% delete(findall(findall(gcf,'Type','axe'),'Type','text'))
%
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%These loops check if the data is gridded, if so it will use the gridded
%gradient matrix calculation. Flagx and Flagy are created. If any of these are 1
%the data is not a uniform grid. The precision is set to 8 (I round) so
%really small grids will suffer here.
%Grid data check X
[ Uniform ] = UniformGridCheck3d( X,Y,Z );

if Uniform==1

    %Calculating spacing
    spacingx=X(1,1,1)-X(1,2,1);
    spacingy=Y(1,1,1)-Y(2,1,1);
    spacingz=Z(1,1,1)-Z(1,1,2);
    

    %If Ux and UY were put in as cols we fix this here
    if iscolumn(Ux)
        Ux=reshape(Ux,size(X));
    end
    if iscolumn(Uy)
        Uy=reshape(Uy,size(X));
    end
    if iscolumn(Uz)
        Uz=reshape(Uz,size(X));
    end

    % Calculating displacement gradients, This is not deformation gradient
    % F, Martel lecture 14 GG303 page 6. See 15 page 3+4
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
    
    progressbar('Calculating Finite strain matrices on non-uniform grid')  
    n=numel(X) ;
    for i=1:n  

        %only getting points closest to first 
        BndA = [X(i),Y(i),Z(i)]; 
        BndB = [X(:),Y(:),Z(:)];
        x2=(BndB(:,1)-BndA(1)).^2; %length x
        y2=(BndB(:,2)-BndA(2)).^2; %length y
        z2=(BndB(:,3)-BndA(3)).^2; %length z
        %total distance from point A
        dist=sqrt(abs(x2)+abs(y2)+abs(z2)); 
        %arranging vec
        BndB = [X(:),Y(:),Z(:),Ux(:),Uy(:),Uz(:)];    
        %sorting on distance, smallest first
        BndB = sortrows([dist,BndB]);        
        %B is now a vector of [distance from A, X vals, Y vals, Zvals, Ux,Uy,Uz]
        %First row can be ignored as it is the same point as A
        %the 20 closest points to the one we are calculating strain on
        closestPoints=BndB(2:21,2:7); 

        %Extracting vars
        C=closestPoints(:,1:3); %[X,Y,Z] coordinates
        d=closestPoints(:,4:6); %[Ux,Uy,Uz]

        % get velocity data and reformat into column vector d
        d = d'; % velocities by row then column
        d = d(:); % make column vector
        % Every 1st row is Ux for the first point and the second is Uy third is Uz 

        % Set up coefficient matrix (Bit cleaner that old method)
        %Zeros
        Zers=zeros(60,1);
        BndA=Zers; BndB=Zers; BndC=Zers;
        Ga=Zers;Gb=Zers;Gc=Zers;
        Gd=Zers;Ge=Zers;Gf=Zers;
        Gg=Zers;Gh=Zers;Gi=Zers;
        %Banded vectors
        BndA(1:3:end)=1; BndB(2:3:end)=1; BndC(3:3:end)=1;
        %Banded vectors of values and 0's with parts of C filling this bands
        Ga(1:3:end)=C(1:1:end,3);
        Gb(1:3:end)=C(1:1:end,2);
        Gc(1:3:end)=C(1:1:end,1);
        Gd(2:3:end)=C(1:1:end,3);
        Ge(2:3:end)=C(1:1:end,2);
        Gf(2:3:end)=C(1:1:end,1);
        Gg(3:3:end)=C(1:1:end,3);
        Gh(3:3:end)=C(1:1:end,2);
        Gi(3:3:end)=C(1:1:end,1);
        %Collate array
        G2=[BndA,BndB,BndC, Ga,Gb,Gc, Gd,Ge,Gf, Gg,Gh,Gi];

        %Solve for model paramaters (tx, ty, tz, Uxz, Uxy, Uxx, Uyz, Uyy, Uyx, Uzz, Uzy, Uzy )
        m=G2\d;
        
        %Filling each row with the calculated unknowns
        Filler(i,:)=m';  

        %If you are getting an error here you have misplaced the file
        %'progressbar.m', just comment the line. 
        progressbar(i/n) % Update figure 
        
    end

    %Extracting parts
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

%Creating a list of displacement tensors, each 3x3 row is a different point. 
n=numel(Ux);
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
F(bad) =nan;
bad = isinf(Ft);
Ft(bad) = nan;

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
%eyx= E(2:N:end,1); %same as exy
eyy= E(2:N:end,2);
eyz= E(2:N:end,3);
%ezx= E(3:N:end,1); %same as ezx
%ezy= E(3:N:end,2); %same as ezy
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
ExxIErr=((exx-exxinf)./exx)*100;
EyyIErr=((eyy-eyyinf)./eyy)*100;
EzzIErr=((ezz-ezzinf)./ezz)*100;
ExyIErr=((exy-exyinf)./exy)*100;
ExzIErr=((exz-exzinf)./exz)*100;
EyxIErr=((eyz-eyzinf)./eyz)*100;
    

end
    
    




 
