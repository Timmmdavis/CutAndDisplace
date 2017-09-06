function [DsTnE1FB,DsTsE1FB,DnTnE1FB,DnTsE1FB,DsTnE1IF,DsTsE1IF,DnTnE1IF,DnTsE1IF,...
Ds_UxE1FB,Ds_UyE1FB,Dn_UxE1FB,Dn_UyE1FB,Ds_UxE1IF,Ds_UyE1IF,Dn_UxE1IF,Dn_UyE1IF,...
DsTnE2FB,DsTsE2FB,DnTnE2FB,DnTsE2FB,DsTnE2IF,DsTsE2IF,DnTnE2IF,DnTsE2IF,...
Ds_UxE2FB,Ds_UyE2FB,Dn_UxE2FB,Dn_UyE2FB,Ds_UxE2IF,Ds_UyE2IF,Dn_UxE2IF,Dn_UyE2IF,...
NUME1,NUME2,nuE1,EE1,nuE2,EE2,...
FreeBoundary1,FreeBoundary2,FdispE1,FdispE2,FreeBoundaries,FixedDisps]...
= InhomogeneousInfluenceBuilder2d...
(x,y,xe,ye,a,Beta,nu,E,halfspace,NormAng,Fdisp)
%Creates and extracts Inhomogeneous influence matrices

%   Copyright 2017, Tim Davis, The University of Aberdeen



    BoundaryFlag=Fdisp; 
    
    %%%%%%%
    %E1 Elastic 1
    %%%%%%%
    Part=1; %Which elastic we want to extract
    [xE1,yE1,xeE1,yeE1,aE1,BetaE1,nuE1,EE1,NormAngE1,NUME1,FdispE1,FB_E1,IF_E1] = ExtractElasticParts2d( Part,BoundaryFlag,x,y,xe,ye,a,Beta,nu,E,NormAng );
    %E1 Elastic 1 creating big influence matrices
    [DsTnE1,DsTsE1,DnTnE1,DnTsE1,Ds_UxE1,Ds_UyE1,Dn_UxE1,Dn_UyE1]...
        = CalculateInfluenceMatrices2d(NUME1,halfspace,xE1,yE1,xeE1,yeE1,aE1,BetaE1,nuE1,EE1,NormAngE1,1 ); %we want disp so passing Fdisp as 1 whatever happens
    clear xE1 yE1 xeE1 yeE1 aE1 BetaE1 NormAngE1   
    %Getting the traction inf matrices for elements on the free boundary and then interface of E1
    [DsTnE1FB,DsTsE1FB,DnTnE1FB,DnTsE1FB] = ExtractData(FB_E1,1,DsTnE1,DsTsE1,DnTnE1,DnTsE1);
    [DsTnE1IF,DsTsE1IF,DnTnE1IF,DnTsE1IF] = ExtractData(IF_E1,1,DsTnE1,DsTsE1,DnTnE1,DnTsE1 );
    clear DsTnE1 DsTsE1 DnTnE1 DnTsE1
    %Getting the displacement inf matrices for elements on the free boundary and then interface of E1    
    [Ds_UxE1FB,Ds_UyE1FB,Dn_UxE1FB,Dn_UyE1FB] = ExtractData(FB_E1,1,Ds_UxE1,Ds_UyE1,Dn_UxE1,Dn_UyE1);    
    [Ds_UxE1IF,Ds_UyE1IF,Dn_UxE1IF,Dn_UyE1IF] = ExtractData(IF_E1,1,Ds_UxE1,Ds_UyE1,Dn_UxE1,Dn_UyE1);     
    clear Ds_UxE1 Ds_UyE1 Dn_UxE1 Dn_UyE1
    
    %%%%%%%
    %E2 Elastic 2
    %%%%%%%
    Part=2; %Which elastic we want to extract
    [xE2,yE2,xeE2,yeE2,aE2,BetaE2,nuE2,EE2,NormAngE2,NUME2,FdispE2,FB_E2,IF_E2]= ExtractElasticParts2d( Part,BoundaryFlag,x,y,xe,ye,a,Beta,nu,E,NormAng );
    %E2 Elastic 2 creating big influence matrices    
    [DsTnE2,DsTsE2,DnTnE2,DnTsE2,Ds_UxE2,Ds_UyE2,Dn_UxE2,Dn_UyE2]...
        = CalculateInfluenceMatrices2d(NUME2,halfspace,xE2,yE2,xeE2,yeE2,aE2,BetaE2,nuE2,EE2,NormAngE2,1); %we want disp so passing Fdisp as 1 whatever happens
    clear halfspace xE2 yE2 xeE2 yeE2 aE2 BetaE2 NormAngE2 
    %Getting the traction inf matrices for elements on the free boundary and then interface of E2    
    [DsTnE2FB,DsTsE2FB,DnTnE2FB,DnTsE2FB] = ExtractData(FB_E2,1,DsTnE2,DsTsE2,DnTnE2,DnTsE2 );
    [DsTnE2IF,DsTsE2IF,DnTnE2IF,DnTsE2IF] = ExtractData(IF_E2,1,DsTnE2,DsTsE2,DnTnE2,DnTsE2 );
    clear DsTnE2 DsTsE2 DnTnE2 DnTsE2
    %Getting the displacement inf matrices for elements on the free boundary and then interface of E2    
    [Ds_UxE2FB,Ds_UyE2FB,Dn_UxE2FB,Dn_UyE2FB] = ExtractData(FB_E2,1,Ds_UxE2,Ds_UyE2,Dn_UxE2,Dn_UyE2 ); 
    [Ds_UxE2IF,Ds_UyE2IF,Dn_UxE2IF,Dn_UyE2IF] = ExtractData(IF_E2,1,Ds_UxE2,Ds_UyE2,Dn_UxE2,Dn_UyE2 ); 
    clear Ds_UxE2 Ds_UyE2 Dn_UxE2 Dn_UyE2        
    

    %Creating some flags for fixing disps
    FreeBoundaries=(BoundaryFlag == 0 | BoundaryFlag == 1 | BoundaryFlag == 4 | BoundaryFlag == 5);%Any free boundaries in either elastic   
    FixedDisps=(BoundaryFlag == 1 | BoundaryFlag == 5);%Any free boundaries in either elastic 
    FixedDisps=FixedDisps(FreeBoundaries); %Fixed displacements on the freeboundary elements
    FreeBoundary2=(BoundaryFlag == 4 | BoundaryFlag == 5);%2nd free boundary
    FreeBoundary2=FreeBoundary2(FreeBoundaries);
    FreeBoundary1=~FreeBoundary2;
    

end


