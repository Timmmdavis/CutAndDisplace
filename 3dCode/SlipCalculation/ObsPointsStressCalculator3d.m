function [TotalStrainStress,StrainStressChange,StrainStressRemote,Sxx,Syy,Szz,Sxy,Sxz,Syz]...
    =ObsPointsStressCalculator3d(mu,lambda,X,Y,Z,Sxx,Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu)
%ObsPointsStressCalculator3d Calculates the strain on points due a fault
%surface with known slip. This is slower than function
%'CalculateStressOnSurroundingPoints' but has the advantage you can store
%these influence matricies if you are running a loop where both the surface
%and obs points geometry stay the same. 

%   Copyright 2017, Tim Davis, The University of Aberdeen

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

%Calculate influence matricies
[DssExx,DssEyy,DssEzz,DssExy,DssExz,DssEyz,DdsExx,DdsEyy,DdsEzz,DdsExy,DdsExz,DdsEyz,...
DnExx,DnEyy,DnEzz,DnExy,DnExz,DnEyz]=CalculateInfluenceMatrices3d_ObsPoints...
(Ss,Ds,Ts,mu,lambda,X,Y,Z,P1,P2,P3,halfspace,nu);

%Accumulate arrays 
Aexx=[DnExx,DssExx,DdsExx];
Aeyy=[DnEyy,DssEyy,DdsEyy];
Aezz=[DnEzz,DssEzz,DdsEzz];
Aexy=[DnExy,DssExy,DdsExy];
Aexz=[DnExz,DssExz,DdsExz];
Aeyz=[DnEyz,DssEyz,DdsEyz];
clear DssExx DssEyy DssEzz DssExy DssExz DssEyz DdsExx DdsEyy DdsEzz DdsExy DdsExz DdsEyz DnExx DnEyy DnEzz DnExy DnExz DnEyz

%Concatenate ready for equation 
A= [Aexx;Aeyy;Aezz;Aexy;Aexz;Aeyz];  
clear Aexx Aeyy Aezz Aexy Aexz Aeyz

%Concatenate vector ready for equation 
B= [Ts;Ss;Ds];  

%Run linear equation system to find result of strains on the obs points
Strains=A*B;

%Now extract variables
NUM = size(X(:),1);
Exx = Strains(1:NUM,:);   
Eyy = Strains((NUM+1):(2*NUM),:);   
Ezz = Strains(((2*NUM)+1):(3*NUM),:);
Exy = Strains(((3*NUM)+1):(4*NUM),:);
Exz = Strains(((4*NUM)+1):(5*NUM),:);
Eyz = Strains(((5*NUM)+1):(6*NUM),:);

%Strain comes out the func, calculating stress from this.  
[Sxx,Syy,Szz,Sxy,Sxz,Syz] = HookesLaw3dStrain2Stress(Exx,...
Eyy,Ezz,Exy,Exz,Eyz,lambda,mu );

%Getting variables ready for output
StrainStressChange=[Exx,Eyy,Ezz,Exy,Exz,Eyz,Sxx,Syy,Szz,Sxy,Sxz,Syz];
TotalStrainStress=StrainStressChange+StrainStressRemote;

end
