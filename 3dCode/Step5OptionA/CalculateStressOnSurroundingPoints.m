function [TotalStrainStress,StrainStressChange,StrainStressRemote,Sxx,Syy,Szz,Sxy,Sxz,Syz]=CalculateStressOnSurroundingPoints(Ss,Ds,Ts,mu,lambda,X,Y,Z,Sxx,Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu)
%   Copyright 2017, Tim Davis, The University of Aberdeen
%CalculateStressOnSurroundingPoints This function imports the fault surfaces with the defined/calculated slip and the observation points. These are run through the stress calculations
%from Mehdi Nikkhoo and Thomas R. Walter published in Geophys. J. Int. (2015) 201, 1119âÿÿ1141. 
%   The calculation output is an array of Xx12 columns where columns1-column6 are exx eyy ezz exy exz eyz and columns7-column12 are Sxx Syy Szz Sxy Sxz Syz. There are 3 outputs: Stress/StrainChange which is the stress
%and strain on the points due to the fault movement. Stress which is the regional stress and strain defined earler by the user and Total Stress which is the regional stress/strain added to the stress change. 
%	A text file is also created called 'StressTensorsTable' that can be used by other software.  


 n = numel(X);
%Placing the defined remote driving stress/strain on arrays that match the size of the XY observation grid points 
 SxxReg=zeros(n,1)+Sxx(1,1);
 SyyReg=zeros(n,1)+Syy(1,1); 
 SzzReg=zeros(n,1)+Szz(1,1); 
 SxyReg=zeros(n,1)+Sxy(1,1); 
 SxzReg=zeros(n,1)+Sxz(1,1); 
 SyzReg=zeros(n,1)+Syz(1,1); 
 
 %Hooke'sLaw to Calcuate strain from the remote stresses that are imported into this function
 [ Exx,Eyy,Ezz,Exy,Exz,Eyz ] = HookesLaw3dStress2Strain( SxxReg,SyyReg,SzzReg,SxyReg,SxzReg,SyzReg,lambda,mu ) ;
  
 %Appending into a big Cvec matrix 
 StrainStressRemote = [Exx,Eyy,Ezz,Exy,Exz,Eyz,SxxReg,SyyReg,SzzReg,SxyReg,SxzReg,SyzReg];
 
%Runs the script to create stress across everypoint in the XYZ grid
%defined. Runs for each triangle and adds the tensors together (superposition). 
 StrainChange = zeros(n,6);
 if halfspace==1
    progressbar('Calculating Change in Stress on Observation XYZ HS') % Create figure and set starting time
    for i = 1:size(P1,1)
    StrainChange = StrainChange + TDstrain_stressHS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Ss(i,:),Ds(i,:),Ts(i,:),mu,lambda);
    progressbar(i/size(P1,1)) % Update figure
    end
 else
    progressbar('Calculating Change in Stress on Observation XYZ FS') % Create figure and set starting time
    for i = 1:size(P1,1)
    StrainChange = StrainChange + TDstrain_stressFS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Ss(i,:),Ds(i,:),Ts(i,:),mu,lambda);
    progressbar(i/size(P1,1)) % Update figure
    end
 end   

%Strain comes out the func, calculating stress from this.  
[SxxChange,SyyChange,SzzChange,SxyChange,SxzChange,SyzChange] = HookesLaw3dStrain2Stress(StrainChange(:,1),...
StrainChange(:,2),StrainChange(:,3),StrainChange(:,4),StrainChange(:,5),StrainChange(:,6),lambda,mu );


StrainStressChange=[StrainChange,SxxChange,SyyChange,SzzChange,SxyChange,SxzChange,SyzChange];

TotalStrainStress=StrainStressChange+StrainStressRemote;
 
% %Spliiting the output tensors into a table with column headers for XYZ and
% %the tensors
% Sxx=reshape(TotalStrainStress(:,7),size(X));
% Syy=reshape(TotalStrainStress(:,8),size(X));
% Szz=reshape(TotalStrainStress(:,9),size(X));
% Sxy=reshape(TotalStrainStress(:,10),size(X));
% Sxz=reshape(TotalStrainStress(:,11),size(X));
% Syz=reshape(TotalStrainStress(:,12),size(X));    %creates coloumns from the TDstress output for each tensor

%Col vecs for export
Sxx = TotalStrainStress(:,7);
Syy = TotalStrainStress(:,8);
Szz = TotalStrainStress(:,9);
Sxy = TotalStrainStress(:,10);
Sxz = TotalStrainStress(:,11);
Syz = TotalStrainStress(:,12);
%StressTensorsTable = table(X2,Y2,Z2,Sxx,Syy,Szz,Sxy,Sxz,Syz); 
%writetable(StressTensorsTable);                  %Writes the table

end 

