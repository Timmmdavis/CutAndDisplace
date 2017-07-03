function [ MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,NUM,Fdisp,FB,IF  ]...
    = ExtractElasticParts3d(Part,BoundaryFlag,MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,Fdisp )
%ExtractElasticParts For Inhomogeneous elastics extracting parts using a
%flag

%   Copyright 2017, Tim Davis, The University of Aberdeen
    %Quick Low Down on BoundaryFlag
    %0=free boundary E1
    %1=fixed bits of the free boundary of E1
    %2=E1-E2 interface, E1 elastic properties, normals point towards E2
    %3=E2-E1 interface, E2 elastic properties, normals point towards E1
    %4=If existed would be free boundary E2, 
    %5=fixed bits of the free boundary of E2
    %Only structured for two elastics at the moment. Introducing E3 would mean there are
    %potentially 6 interfaces. E1-E2 E1-E3, E2-E1 E2=E3, E3-E1 E3-E2
    
    if Part==1
    Flag=BoundaryFlag<=2;% Flag for E1   
    BoundaryFlag=BoundaryFlag(Flag);%E1 parts
    
    FB=(BoundaryFlag == 0 | BoundaryFlag == 1);%size of free boundary on E1 (including fixed bits )
    Fdisp=BoundaryFlag==1;%Fdisp on E1        
    IF=BoundaryFlag==2;%size of interface on E1

    
    elseif Part==2
    Flag=BoundaryFlag>=3;% Flag for E2    
    BoundaryFlag=BoundaryFlag(Flag);%E2 parts
    
    IF=BoundaryFlag==3;%size of interface boundary on E2    
    FB=(BoundaryFlag == 4 | BoundaryFlag == 5);%size of free boundary on E2 (including fixed bits )
    Fdisp=BoundaryFlag==5;%Fdisp on E2        

    
    end
    
    NUM = sum(Flag);
    nu=nu(Part);
    mu=mu(Part);
    lambda=lambda(Part);

    
    P1=P1([Flag,Flag,Flag]);     
    P2=P2([Flag,Flag,Flag]);  
    P3=P3([Flag,Flag,Flag]);      
    MidPoint=MidPoint([Flag,Flag,Flag]);
    FaceNormalVector=FaceNormalVector([Flag,Flag,Flag]);
    P1=reshape(P1,[],3);
    P2=reshape(P2,[],3);
    P3=reshape(P3,[],3);    
    MidPoint=reshape(MidPoint,[],3);
    FaceNormalVector=reshape(FaceNormalVector,[],3);

end

