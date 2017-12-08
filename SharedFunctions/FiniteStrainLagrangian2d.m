function [exx,eyy,exy,ExxIErr,EyyIErr,ExyIErr]=FiniteStrainLagrangian2d(Ux,Uy,X,Y)
% FiniteStrainLagrangian2d: Calculate the Lagrangian strain from input grid 
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
% [exx,eyy,exy,ExxIErr,EyyIErr,ExyIErr]=FiniteStrainLagrangian2d(Ux,Uy,X,Y)
%
% usage #2: If the infinitesimal strain errors dont concern you. 
% [exx,eyy,exy]=FiniteStrainLagrangian2d(Ux,Uy,X,Y)
%
% Arguments: (input)
% Ux,Uy             - Displacement components in X and Y at points defined
%                    in the vectors X and Y.
%
% X,Y               - The location of the points in 2D space.
%
% Arguments: (output)
% exx,eyy,exy       - The 2-D Lagrangian finite strain tensors,
%                    if you want infinitesimal ones these are computed also
%                    and could be exported.
%
% ExxIErr,EyyIErr
% ExyIErr           - Error in percent if we assume infinitesimal strain. 
%
%
% Example usage 1: Gridded data:
% 
%  [X,Y]=meshgrid(-2:0.1:2,-2:0.1:2);
%  [~,r] = cart2pol(X,Y);
%  Ux=(r(:).^2).*X(:);Uy=(r(:).^2).*Y(:);
%  [exx,eyy,exy]=FiniteStrainLagrangian2d(Ux,Uy,X,Y);
%  Dilatation=(exx+eyy)/2;
%  Scl=0.1;
%  cmap='default';
%  DrawDeformedGrid2d( X,Y,Ux,Uy,cmap,Dilatation,'Scale',0.1 );
% 
% Example usage 2: Non gridded data:
%
% Density=5000;
% n=4;
% xmv=-2; ymv=-2;zmv=0;
% X=(rand(1,Density)*n)+xmv;
% Y=(rand(1,Density)*n)+ymv;
% [~,r] = cart2pol(X,Y);
% Ux=(r.^2).*X;Uy=(r.^2).*Y;
% [exx,eyy,exy]=FiniteStrainLagrangian2d(Ux,Uy,X,Y);
% Dilatation=(exx+eyy)/2;
% Scl=0.1;
% cmap='default';
% DrawScatterPlots2d( (X+(Ux*Scl)),(Y+(Uy*Scl)),cmap, Dilatation )
% caxis([min(Dilatation(:)) max(Dilatation(:))])
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%These loops check if the data is gridded, if so it will use the gridded
%gradient matrix calculation. Flagx and Flagy are created. If any of these are 1
%the data is not a uniform grid. The precision is set to 8 (I round) so
%really small grids will suffer here.
%Grid data check X
[ Uniform ] = UniformGridCheck2d( X,Y );

if Uniform==1

    %Calculating spacing
    spacingx=X(1,1)-X(1,2);
    spacingy=Y(1,1)-Y(2,1);
    
    %If Ux and UY were put in as cols we fix this here
    if iscolumn(Ux)
        Ux=reshape(Ux,size(X));
    end
    if iscolumn(Uy)
        Uy=reshape(Uy,size(X));
    end

    % Calculating displacement gradients, This is not deformation gradient
    % F, Martel lecture 14 GG303 page 6. See 15 page 3+4
    [Uxx,Uxy] = gradient(Ux,abs(spacingx)); %changed from D.Kroons to match Steve Martels lecture notes
    [Uyx,Uyy] = gradient(Uy,abs(spacingy));
    
else %could use interpolation instead of this method but that can give dodgy edge gradients. Neither are perfect
    
    X=X(:);
    Y=Y(:);
    Ux=Ux(:);
    Uy=Uy(:);
    Filler=zeros(numel(X),6);

    progressbar('Calculating Finite strain matrices on non-uniform grid')  
    n=numel(X) ;
    for i=1:n
       
        %only getting points closest to first 
        BndA = [X(i),Y(i)]; 
        BndB = [X(:),Y(:)];
        %total distance from point A
        dist=bsxfun(@hypot,BndB(:,1)-BndA(1),BndB(:,2)-BndA(2)); 
        %arranging vec
        BndB = [X(:),Y(:),Ux(:),Uy(:)];   
        %sorting on distance, smallest first
        BndB = sortrows([dist,BndB]);         
        %B is now a vector of [distance from A, X vals, Y vals,Ux,Uy]
        %First row can be ignored as it is the same point as A
        %the 20 closest points to the one we are calculating strain on
        closestPoints=BndB(2:21,2:5); 

        %Extracting vars
        C=closestPoints(:,1:2); %[X,Y] coordinates
        d=closestPoints(:,3:4); %[Ux,Uy]

        % get velocity data and reformat into column vector d
        d = d'; % velocities by row then column
        d = d(:); % make column vector
        % Every 1st row is Ux for the first point and the second is Uy 

        % Set up coefficient matrix (Bit cleaner that old method)
        %Zeros
        Zers=zeros(40,1);
        BndA=Zers; BndB=Zers;
        Ga=Zers;Gb=Zers;
        Gc=Zers;Gd=Zers;
        %Banded vectors
        BndA(1:2:end)=1; BndB(2:2:end)=1;
        %Banded vectors of values and 0's with parts of C filling this bands
        Ga(1:2:end)=C(1:1:end,2);
        Gb(1:2:end)=C(1:1:end,1);
        Gc(2:2:end)=C(1:1:end,2);
        Gd(2:2:end)=C(1:1:end,1);
        %Collate array        
        G2=[BndA,BndB, Ga,Gb, Gc,Gd];

        %Solve for model paramaters (tx, ty, Uxy, Uxx, Uyy, Uyx)
        m=G2\d;

        %Filling each row with the calculated unknowns
        Filler(i,:)=m'; 

        %If you are getting an error here you have misplaced the file
        %'progressbar.m', just comment the line. 
        progressbar(i/n) % Update figure 
    
    end

    %Extracting parts
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
F(bad) = nan;
bad = isinf(Ft);
Ft(bad) = nan;

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
%eyx= E(2:N:end,1); %same as exy
eyy= E(2:N:end,2);

%Using calculated infinitesimal strain to calcluate error associated
%with infinitesimal calcluation compared to finite. 
%See D Pollard  2005 Chapter 5.5.3
exxinf= e(1:N:end,1);
exyinf= e(1:N:end,2);
eyyinf= e(2:N:end,2);
ExxIErr=abs(((exx-exxinf)./exx)*100);
EyyIErr=abs(((eyy-eyyinf)./eyy)*100);
ExyIErr=abs(((exy-exyinf)./exy)*100);


end
 
