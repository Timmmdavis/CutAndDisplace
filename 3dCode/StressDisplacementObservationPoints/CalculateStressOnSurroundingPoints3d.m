function [StressTTotal,StrainTTotal,StressTChg,StrainTChg,StressTReg,StrainTReg]...
    =CalculateStressOnSurroundingPoints3d(Dss,Dds,Dn,mu,...
    lambda,X,Y,Z,Sxx,Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu)
% CalculateStressOnSurroundingPoints: calculates the displacement by
%				superposition.
%               This runs the calculated/defined slip for each element and
%               works out the induced disp in the surrounding points. The
%               disps for each element are summed to a the total disp
%               for each point within the loop.
%
%               Commented lines at the base of the script allow writing of
%               a table of the stress tensors. 
%
% usage #1:
%[StressTTotal,StrainTTotal,StressTChg,StrainTChg,StressTReg,StrainTReg]...
% =CalculateStressOnSurroundingPoints3d(Dss,Dds,Dn,mu,...
%    lambda,X,Y,Z,Sxx,Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu)
%
% Arguments: (input)
%    X,Y & Z  - The observation points
%
% P1,P2,P3    - n*3 Column vectors where each 'P' represents the
%               different corner points of one of the triangles (XYZ).
%
%  Dss,Dds,Dn - Vectors that describe how much the elements displace in the
%               normal (Dn) and strike slip (Dss) and dipslip (Dds)
%               directions on the elements.
%
% Sxx,Syy,Szz
% Sxy,Sxz,Syz- Remote stress if it existed on import (useful if the
%              user defined strain boundary conditions. 
%
%       mu    - Shear modulus.
%
%     lambda  - Lame's constant.
%
%       nu    - The Poisson's ratio.
%
%
%  halfspace  - Defines if we work out the coefficientsin a half or whole
%              space
%
% Arguments: (output)
% StressTTotal    - StressTensorTotal. From the stress calculation 6*n 
%                   vector,[Sxx,Syy,Szz,Sxy,Sxz,Syz]. Remote stress and
%                   stress change induced by movement of dislocations.
%
% StrainTTotal    - StrainTensorTotal From the strain calculation 6*n 
%                   vector, [Exx,Eyy,Ezz,Exy,Exz,Eyz]. Remote strain and
%                   strain change induced by movement of dislocations.
%
% StressTChg      - StressTensorChange. From the stress calculation 6*n 
%                   vector, [Sxx,Syy,Szz,Sxy,Sxz,Syz]. Stress change
%                   induced by movement of dislocations.
%
% StrainTChg      - StrainTensorChange From the strain calculation 6*n 
%                   vector, [Exx,Eyy,Ezz,Exy,Exz,Eyz]. Strain change
%                   induced by movement of dislocations.
%
% StressTReg      - StressTensorRegionalStress. From the stress calculation
%                   6*n vector, [Sxx,Syy,Szz,Sxy,Sxz,Syz]. Remote stress
%                   (far field tectonic, not gravity, this would need to be
%                   added independently).
%
% StrainTReg      - StressTensorRegionalStrain From the strain calculation 
%                   6*n vector, [Exx,Eyy,Ezz,Exy,Exz,Eyz]. Remote strain
%                   (far field tectonic)
%
% Example usage:
% 
% N/A
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Checking elastic constants match
ElasticConstantsCheck( mu,lambda,nu );

n = numel(X);

%Placing the defined remote driving stress/strain on arrays that match the size of the XY observation grid points 
SxxReg=zeros(n,1)+Sxx(1,1);
SyyReg=zeros(n,1)+Syy(1,1); 
SzzReg=zeros(n,1)+Szz(1,1); 
SxyReg=zeros(n,1)+Sxy(1,1); 
SxzReg=zeros(n,1)+Sxz(1,1); 
SyzReg=zeros(n,1)+Syz(1,1); 
%Hooke'sLaw to Calcuate strain from the remote stresses that are imported into this function
[ ExxReg,EyyReg,EzzReg,ExyReg,ExzReg,EyzReg ] = HookesLaw3dStress2Strain( SxxReg,SyyReg,SzzReg,SxyReg,SxzReg,SyzReg,lambda,mu ) ;

%Appending into a big Cvec matrix 
StrainTReg=[ExxReg,EyyReg,EzzReg,ExyReg,ExzReg,EyzReg];
StressTReg=[SxxReg,SyyReg,SzzReg,SxyReg,SxzReg,SyzReg];

%Runs the script to create stress across everypoint in the XYZ grid
%defined. Runs for each triangle and adds the tensors together (superposition). 
StrainTChg = zeros(n,6);
if halfspace==1
    progressbar('Calculating Change in Stress on Observation XYZ HS') % Create figure and set starting time
    for i = 1:size(P1,1)
        StrainTChg = StrainTChg + TDstrainHS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Dss(i,:),Dds(i,:),Dn(i,:),mu,lambda);
        progressbar(i/size(P1,1)) % Update figure
    end
    else
    progressbar('Calculating Change in Stress on Observation XYZ FS') % Create figure and set starting time
    for i = 1:size(P1,1)
        StrainTChg = StrainTChg + TDstrainFS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Dss(i,:),Dds(i,:),Dn(i,:),mu,lambda);
        progressbar(i/size(P1,1)) % Update figure
    end
end   

%Strain comes out the func, calculating stress from this.  
[SxxChange,SyyChange,SzzChange,SxyChange,SxzChange,SyzChange] = HookesLaw3dStrain2Stress(StrainTChg(:,1),...
StrainTChg(:,2),StrainTChg(:,3),StrainTChg(:,4),StrainTChg(:,5),StrainTChg(:,6),lambda,mu );

StressTChg=[SxxChange,SyyChange,SzzChange,SxyChange,SxzChange,SyzChange];

StressTTotal=StressTChg+StressTReg;
StrainTTotal=StrainTChg+StrainTReg;
 
%Col vecs for export
% Sxx = StressTTotal(:,1);
% Syy = StressTTotal(:,2);
% Szz = StressTTotal(:,3);
% Sxy = StressTTotal(:,4);
% Sxz = StressTTotal(:,5);
% Syz = StressTTotal(:,6);
%StressTensorsTable = table(X2,Y2,Z2,Sxx,Syy,Szz,Sxy,Sxz,Syz); 
%writetable(StressTensorsTable);                  %Writes the table

end 

