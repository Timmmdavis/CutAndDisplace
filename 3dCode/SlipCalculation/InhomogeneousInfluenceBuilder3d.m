function [DnTnE1FB,DnTssE1FB,DnTdsE1FB,DssTnE1FB,DssTssE1FB,DssTdsE1FB,DdsTnE1FB,DdsTssE1FB,DdsTdsE1FB,...
DnTnE1IF,DnTssE1IF,DnTdsE1IF,DssTnE1IF,DssTssE1IF,DssTdsE1IF,DdsTnE1IF,DdsTssE1IF,DdsTdsE1IF,...
Dn_dxE1FB,Dn_dyE1FB,Dn_dzE1FB,Dss_dxE1FB,Dss_dyE1FB,Dss_dzE1FB,Dds_dxE1FB,Dds_dyE1FB,Dds_dzE1FB,...
Dn_dxE1IF,Dn_dyE1IF,Dn_dzE1IF,Dss_dxE1IF,Dss_dyE1IF,Dss_dzE1IF,Dds_dxE1IF,Dds_dyE1IF,Dds_dzE1IF,...
DnTnE2FB,DnTssE2FB,DnTdsE2FB,DssTnE2FB,DssTssE2FB,DssTdsE2FB,DdsTnE2FB,DdsTssE2FB,DdsTdsE2FB,...
DnTnE2IF,DnTssE2IF,DnTdsE2IF,DssTnE2IF,DssTssE2IF,DssTdsE2IF,DdsTnE2IF,DdsTssE2IF,DdsTdsE2IF,...
Dn_dxE2FB,Dn_dyE2FB,Dn_dzE2FB,Dss_dxE2FB,Dss_dyE2FB,Dss_dzE2FB,Dds_dxE2FB,Dds_dyE2FB,Dds_dzE2FB,...
Dn_dxE2IF,Dn_dyE2IF,Dn_dzE2IF,Dss_dxE2IF,Dss_dyE2IF,Dss_dzE2IF,Dds_dxE2IF,Dds_dyE2IF,Dds_dzE2IF,...
NUME1,NUME2,nuE1,muE1,nuE2,muE2,...
FreeBoundary1,FreeBoundary2,FdispE1,FdispE2,FreeBoundaries,FixedDisps]...
= InhomogeneousInfluenceBuilder3d...
(MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,halfspace,nu,Fdisp)
%Creates and extracts Inhomogeneous influence matrices

%   Copyright 2017, Tim Davis, The University of Aberdeen
    BoundaryFlag=Fdisp; 
    
    %%%%%%%
    %E1 Elastic 1
    %%%%%%%
    Part=1; %Which elastic we want to extract
    [MidPointE1,P1E1,P2E1,P3E1,muE1,lambdaE1,FaceNormalVectorE1,nuE1,NUME1,FdispE1,FB_E1,IF_E1]...
    = ExtractElasticParts3d( Part,BoundaryFlag,MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,Fdisp );
    
    [ DnTnE1,DnTssE1,DnTdsE1,DssTnE1,DssTssE1,DssTdsE1,DdsTnE1,DdsTssE1,DdsTdsE1,Dn_dxE1,Dn_dyE1,...
    Dn_dzE1,Dss_dxE1,Dss_dyE1,Dss_dzE1,Dds_dxE1,Dds_dyE1,Dds_dzE1,NUME1,StrikeSlipCosineE1,DipSlipCosineE1]...
    = CalculateInfluenceMatrices3d(MidPointE1,P1E1,P2E1,P3E1,muE1,lambdaE1,FaceNormalVectorE1,halfspace,nuE1,1);
    clear MidPointE1 P1E1 P2E1 P3E1  lambdaE1   

    %Getting the traction inf matrices for elements on the free boundary and then interface of E1    
    [DnTnE1FB,DnTssE1FB,DnTdsE1FB,DssTnE1FB,DssTssE1FB,DssTdsE1FB,DdsTnE1FB,DdsTssE1FB,DdsTdsE1FB]...
    = ExtractData(FB_E1,1,DnTnE1,DnTssE1,DnTdsE1,DssTnE1,DssTssE1,DssTdsE1,DdsTnE1,DdsTssE1,DdsTdsE1 );
    [DnTnE1IF,DnTssE1IF,DnTdsE1IF,DssTnE1IF,DssTssE1IF,DssTdsE1IF,DdsTnE1IF,DdsTssE1IF,DdsTdsE1IF]...
    = ExtractData(IF_E1,1,DnTnE1,DnTssE1,DnTdsE1,DssTnE1,DssTssE1,DssTdsE1,DdsTnE1,DdsTssE1,DdsTdsE1 );
    clear  DnTnE1 DnTssE1 DnTdsE1 DssTnE1 DssTssE1 DssTdsE1 DdsTnE1 DdsTssE1 DdsTdsE1
    %Getting the displacement inf matrices for elements on the free boundary and then interface of E1 
    [Dn_dxE1FB,Dn_dyE1FB,Dn_dzE1FB,Dss_dxE1FB,Dss_dyE1FB,Dss_dzE1FB,Dds_dxE1FB,Dds_dyE1FB,Dds_dzE1FB]...
    = ExtractData(FB_E1,1 ,Dn_dxE1,Dn_dyE1,Dn_dzE1,Dss_dxE1,Dss_dyE1,Dss_dzE1,Dds_dxE1,Dds_dyE1,Dds_dzE1);
    [Dn_dxE1IF,Dn_dyE1IF,Dn_dzE1IF,Dss_dxE1IF,Dss_dyE1IF,Dss_dzE1IF,Dds_dxE1IF,Dds_dyE1IF,Dds_dzE1IF]...
    = ExtractData(IF_E1,1,Dn_dxE1,Dn_dyE1,Dn_dzE1,Dss_dxE1,Dss_dyE1,Dss_dzE1,Dds_dxE1,Dds_dyE1,Dds_dzE1);
    clear  Dn_dxE1 Dn_dyE1 Dn_dzE1 Dss_dxE1 Dss_dyE1 Dss_dzE1 Dds_dxE1 Dds_dyE1 Dds_dzE1

    
    %%%%%%%
    %E2 Elastic 2
    %%%%%%%
    Part=2; %Which elastic we want to extract
    [MidPointE2,P1E2,P2E2,P3E2,muE2,lambdaE2,FaceNormalVectorE2,nuE2,NUME2,FdispE2,FB_E2,IF_E2]...
    = ExtractElasticParts3d( Part,BoundaryFlag,MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,Fdisp );
    
    [ DnTnE2,DnTssE2,DnTdsE2,DssTnE2,DssTssE2,DssTdsE2,DdsTnE2,DdsTssE2,DdsTdsE2,Dn_dxE2,Dn_dyE2,...
    Dn_dzE2,Dss_dxE2,Dss_dyE2,Dss_dzE2,Dds_dxE2,Dds_dyE2,Dds_dzE2,NUME2,StrikeSlipCosineE2,DipSlipCosineE2]...
    = CalculateInfluenceMatrices3d(MidPointE2,P1E2,P2E2,P3E2,muE2,lambdaE2,FaceNormalVectorE2,halfspace,nuE2,1);

    clear MidPointE2 P1E2 P2E2 P3E2  lambdaE2    

    
            %Getting the traction inf matrices for elements on the free boundary and then interface of E2    
    [DnTnE2FB,DnTssE2FB,DnTdsE2FB,DssTnE2FB,DssTssE2FB,DssTdsE2FB,DdsTnE2FB,DdsTssE2FB,DdsTdsE2FB]...
    = ExtractData(FB_E2,1,DnTnE2,DnTssE2,DnTdsE2,DssTnE2,DssTssE2,DssTdsE2,DdsTnE2,DdsTssE2,DdsTdsE2 );
    [DnTnE2IF,DnTssE2IF,DnTdsE2IF,DssTnE2IF,DssTssE2IF,DssTdsE2IF,DdsTnE2IF,DdsTssE2IF,DdsTdsE2IF]...
    = ExtractData(IF_E2,1,DnTnE2,DnTssE2,DnTdsE2,DssTnE2,DssTssE2,DssTdsE2,DdsTnE2,DdsTssE2,DdsTdsE2 );
    clear  DnTnE2 DnTssE2 DnTdsE2 DssTnE2 DssTssE2 DssTdsE2 DdsTnE2 DdsTssE2 DdsTdsE2
    %Getting the displacement inf matrices for elements on the free boundary and then interface of E2 
    [Dn_dxE2FB,Dn_dyE2FB,Dn_dzE2FB,Dss_dxE2FB,Dss_dyE2FB,Dss_dzE2FB,Dds_dxE2FB,Dds_dyE2FB,Dds_dzE2FB]...
    = ExtractData(FB_E2,1,Dn_dxE2,Dn_dyE2,Dn_dzE2,Dss_dxE2,Dss_dyE2,Dss_dzE2,Dds_dxE2,Dds_dyE2,Dds_dzE2 );
    [Dn_dxE2IF,Dn_dyE2IF,Dn_dzE2IF,Dss_dxE2IF,Dss_dyE2IF,Dss_dzE2IF,Dds_dxE2IF,Dds_dyE2IF,Dds_dzE2IF]...
    = ExtractData(IF_E2,1,Dn_dxE2,Dn_dyE2,Dn_dzE2,Dss_dxE2,Dss_dyE2,Dss_dzE2,Dds_dxE2,Dds_dyE2,Dds_dzE2);
    clear  Dn_dxE2 Dn_dyE2 Dn_dzE2 Dss_dxE2 Dss_dyE2 Dss_dzE2 Dds_dxE2 Dds_dyE2 Dds_dzE2
    
    
    %Creating some flags for fixing disps
    FreeBoundaries=(BoundaryFlag == 0 | BoundaryFlag == 1 | BoundaryFlag == 4 | BoundaryFlag == 5);%Any free boundaries in either elastic   
    FixedDisps=(BoundaryFlag == 1 | BoundaryFlag == 5);%Any free boundaries in either elastic 
    FixedDisps=FixedDisps(FreeBoundaries); %Fixed displacements on the freeboundary elements
    FreeBoundary2=(BoundaryFlag == 4 | BoundaryFlag == 5);%2nd free boundary
    FreeBoundary2=FreeBoundary2(FreeBoundaries);
    FreeBoundary1=~FreeBoundary2;
    
    

    


end


