function [Dss,Dds,Dn,Pxx,Pyy,Pzz,Pxy,Pxz,Pyz]=SlipCalculator3d(MidPoint,...
Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,TnIn,TssIn,TdsIn,mu,lambda,nu,P1,P2,P3,halfspace...
,FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points,cmap)
% SlipCalculator2d: Calculates the 3 slip components on the
%               imported fault surface from  boundary condition defined by the
%               user.
%
% usage #1:
%
% [Dss,Dds,Dn,Sxx,Syy,Szz,Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
% Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option);
%
% Arguments: (input)
%    MidPoint - The element midpoints (triangle). Incentre.  
%
% Pxx,Pyy,Pzz
% Pxy,Pxz,Pyz - Remote stress defined by user.
%
% TnIn,TssIn
% TdsIn       - Tractions on the surface defined by the user 
%              (Normal,Strikeslip and Dipslip) .
%
%       mu    - Shear modulus.
%
%     lambda  - Lame's constant.
%
%       nu    - The Poisson's ratio.
%
% P1,P2,P3    - n*3 Column vectors where each 'P' represents the
%               different corner points of one of the triangles (XYZ).
%
%  halfspace  - Defines if we work out the coefficientsin a half or whole
%              space
%
% FaceNormalVector - The direction cosines, CosAx (Nx), CosAy and CosAz in 
%                   a list for each element. 
%
%    Strain   - Flag to say if the input stresses are actually a strain
%              tensor. If so we convert these to stresses before
%              continuing.
%
% Fdisp       - Flag telling the user if any elements are going to be fixed
%              (if this is the case displacement influence matricies are
%              required). In the case of inhomogeneous materials this is
%              more complex, see function
%              'InhomogeneousInfluenceBuilder3d.m'.
%
%   Mu        - Coefficient of friction at each element.
%
%   Sf        - Sliding friction (Cohesion) at each element.
%
%   Option    - which option the user is using:
%           'B'=under remote stress, fault can interpenetrate.
%           'C'=friction solver as in Ritz 2012.
%           'D'=remote stress, opening components not calculated. Fault satisfies
%             boundary conditions by shearing only.
%           'E'=tractions defined on body, I.e. pressure in a magma
%           chamber.
%           'F'=Inhomogeneous.
%
% Arguments: (output)
%  Dss,Dds,Dn - Vectors that describe how much the elements displace in the
%               normal (Dn) and strike slip (Dss) and dipslip (Dds)
%               directions on the elements.
%
% Pxx,Pyy,Pzz
% Pxy,Pxz,Pyz- Remote stress if it existed on import (useful if the
%              user defined strain boundary conditions. 
%
% Example usage:
% %Assuming the elastic constants and surface are defined:
%  halfspace=0;
%  [Sxx,Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,Mu,Sf,strain,SecondSurface ]...
%  = CreateBlankVars;
%  strain=0; 
%  Sxx = 0; 				
%  Syy = 0; 
%  Szz = 0; 
%  Sxy = 1; 					
%  Sxz = 0; 					
%  Syz = 0;                
%  Option='B'; 
%  [Dss,Dds,Dn,Sxx,Syy,Szz,Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
%   Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option);
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Number of points:
NUM=numel(MidPoint(:,1));

%Check there are no duplicate els. 
DuplicateElementCheck(MidPoint,Fdisp)

%Check input constants match:
ElasticConstantsCheck( mu(1),lambda(1),nu(1) );
if numel(nu)>1 %Inhomo, 2 elastics
    ElasticConstantsCheck( mu(2),lambda(2),nu(2) );
end


% If Defined slip is as a uniform remote strain, checking options that have
% a uniform 'stress' defined
if Option=='B' || Option=='C' || Option=='D' || Option=='F'

    %for the imported stress making sure these are col vectors
    [ Pxx,Pyy,Pzz,Pxy,Pxz,Pyz ] = RowVecToCol( Pxx,Pyy,Pzz,Pxy,Pxz,Pyz );
    
    if strain==1
        %Converting the input strain boundary condition on every faults midpoint to stress
        %Hooke's law, Equation 7.131 and 7.132 in David Pollards Book.  
        [ Pxx,Pyy,Pzz,Pxy,Pxz,Pyz ] = HookesLaw3dStrain2Stress( Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,lambda,mu );
    end  

end

%Function to check the boundary stress tensors related to Z(y in 2d) components are 0 at the
%halfspace surface
if halfspace==1
    HalfSpaceBoundaryConditionsCheck(MidPoint(:,3),Pzz,Pyz,Pxz)
end

%Checking we can overcome the frictional resistance on the fault. No
%point waiting for this otherwise
if Option=='C'     
    
    [ Tn,Tds,Tss ] = CalculateNormalAndShearTractions3d( FaceNormalVector,Pxx,Pyy,Pzz,Pxy,Pxz,Pyz );
    %Adding tractions imported into function if these also exist. 
    if ~isempty(TnIn) && ~isempty(TdsIn) && ~isempty(TssIn)    
        Tn=Tn+TnIn;
        Tds=Tds+TdsIn;
        Tss=Tss+TssIn;
    end
    FricRes=Tn.*Mu; %Negative is compression
    if all(-abs(Tss)>FricRes) && all(-abs(Tds)>FricRes)
        cmap2 = colormap_cpt('Ccool-warm2');
        PlotSlipDistribution3d(Triangles,Points,cmap2,FricRes,Tss,Tds)
        disp('Your frictional fault will not slip, quit (ctrl+c)')
        pause        
    end
    
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating influence matrices using a function. Big square matrices of
%how much every element slipping by 1 effects traction on every other element
%If displacement if fixed on any elements we also disp influence matrices. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if Option~='F'  %any option but inhomo   
    
    [StressInf,DispInf] = CalculateInfluenceMatrices3d(MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,halfspace,nu,Fdisp);
    
else %Inhomo calc
    
    [StressInfE1FB,StressInfE1IF,DispInfE1FB,DispInfE1IF,StressInfE2FB,StressInfE2IF,DispInfE2FB,DispInfE2IF,...
     NUME1,NUME2,FreeBoundary1,FreeBoundary2,FdispE1,FdispE2,FreeBoundaries,FixedDisps]...
     = InhomogeneousInfluenceBuilder3d...
     (MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,halfspace,nu,Fdisp);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Traction influence matrices should now exist. 
%Now the defined pressures, Dss,Dds and Dn are reshaped. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Option=='B' || Option=='C' || Option=='D' || Option=='F' 
    
    if Option=='F'   

        FaceNormalVector=FaceNormalVector([FreeBoundaries,FreeBoundaries,FreeBoundaries]); 
        FaceNormalVector=reshape(FaceNormalVector,[],3); 
        
    end     
    
    [ Tn,Tds,Tss ] = CalculateNormalAndShearTractions3d( FaceNormalVector,Pxx,Pyy,Pzz,Pxy,Pxz,Pyz );
    
    %Adding tractions imported into function if these also exist. 
    if ~isempty(TnIn) && ~isempty(TdsIn) && ~isempty(TssIn)    
		disp('Boundarys have both tractions and remote stresses defined.')
        Tn=Tn+TnIn;
        Tds=Tds+TdsIn;
        Tss=Tss+TssIn;
    end
    
end

    %for option E traction has been already defined
if Option=='E' 
    
    %If Tractions do not exist yet (i.e. option E) then we create blank ones
    %to append the user defined ones too. 
    Tn=zeros(NUM,1);
    Tss=zeros(NUM,1);
    Tds=zeros(NUM,1);
    
    %creates a column of the stress on every triangles centre point with 
    %the pressure direction defined by the slip type in the name. 
    Tn      = Tn+TnIn;    
    Tss     = Tss+TssIn;
    Tds     = Tds+TdsIn;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating large arrays from the influence matrices and the traction matrices.
%performing the system of linear equations to find the slip due to the boundary conditions. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Quick loop to force any fixed triangles to stick
if Option~='F'  
    
    if any(Fdisp)==1
        %Now exchanging stress coeff rows for displacement ones
        [StressInf,NUM,Tn,Tss,Tds]=FixingDisp_InfRowSwap3d(StressInf,DispInf,NUM,Tn,Tss,Tds,Fdisp);
        clear DispInf
    end
    
else %We need to do this for both elastics  
    
    if any(FdispE1)==1    
        [ StressInfE1FB,NUME1,~,~,~]= FixingDisp_InfRowSwap3d( StressInfE1FB,DispInfE1FB,NUME1,[],[],[],FdispE1);       
        %Now removing any fixed disp element cols from the interface elements  
        [StressInfE1IF] = ExtractData(FdispE1,2,StressInfE1IF,'Columns',1 );
        [DispInfE1IF]  = ExtractData(FdispE1,2,DispInfE1IF,'Columns',1 );
    end
    
    %E2    
    if any(FdispE2)==1
        %Need to Flip So the fixed disps are the top part    
        [StressInfE2FB,DispInfE2FB,FdispE2]=FlipColVecs(StressInfE2FB,DispInfE2FB,FdispE2);
        
        [ StressInfE2FB,NUME1,~,~,~ ]= FixingDisp_InfRowSwap3d( StressInfE2FB,DispInfE2FB,NUME2,[],[],[],FdispE2);  
        
        %Flipping back right way up
        [StressInfE2FB,DispInfE2FB,FdispE2]=FlipColVecs(StressInfE2FB,DispInfE2FB,FdispE2);
        %Now removing any fixed disp element cols from the interface elements  
        [StressInfE2IF] = ExtractData(FdispE2,2,StressInfE2IF,'Columns',1 );
        [DispInfE2IF]  = ExtractData(FdispE2,2,DispInfE2IF,'Columns',1 );      
    end
    
    FD=logical(FixedDisps);    
    Tn(FD) = 0;    
    Tss(FD) = 0;
    Tds(FD) = 0;
    
end       
    
if Option=='B' || Option=='E' || Option=='C'
        
    Atn  = [-StressInf.DnTn,  -StressInf.DssTn, -StressInf.DdsTn ];     StressInf = rmfield(StressInf,{'DnTn','DssTn','DdsTn'});
    Atss = [-StressInf.DnTss, -StressInf.DssTss,-StressInf.DdsTss];     StressInf = rmfield(StressInf,{'DnTss','DssTss','DdsTss'});
    Atds = [-StressInf.DnTds, -StressInf.DssTds,-StressInf.DdsTds];     
    clear StressInf
    A= [Atn;Atss;Atds];  %Concatenate ready for equation 
    clear Atn Atss Atds

elseif Option=='D'             
                                   
    Atn =  [-StressInf.DssTn, -StressInf.DdsTn ];                       StressInf = rmfield(StressInf,{'DssTn','DdsTn'});
    Atss = [-StressInf.DssTss,-StressInf.DdsTss];                       StressInf = rmfield(StressInf,{'DssTss','DdsTss'});
    Atds = [-StressInf.DssTds,-StressInf.DdsTds];
    clear StressInf
    A= [Atn;Atss;Atds];  %Concatenate ready for equation 
    clear Atn Atss Atds
    
elseif Option=='F'         
    %Note this is not currently very efficient way of doing the memory allocation. Better to build each row of 'A' separately.  
    zerinfA=zeros(size(StressInfE1FB.DnTn,1),size(DispInfE2IF.DnUz,2)); %Sz=1=(Rows down),2=(Cols across)
    zerinfB=zeros(size(StressInfE2FB.DnTn,1),size(DispInfE1IF.DnUz,2)); 
    
    
    A=[[StressInfE1IF.DnTn,   StressInfE1IF.DssTn,  StressInfE1IF.DdsTn,  -StressInfE2IF.DnTn,  -StressInfE2IF.DssTn,  -StressInfE2IF.DdsTn ];
       [StressInfE1IF.DnTss,  StressInfE1IF.DssTss, StressInfE1IF.DdsTss, -StressInfE2IF.DnTss, -StressInfE2IF.DssTss, -StressInfE2IF.DdsTss];
       [StressInfE1IF.DnTds,  StressInfE1IF.DssTds, StressInfE1IF.DdsTds, -StressInfE2IF.DnTds, -StressInfE2IF.DssTds, -StressInfE2IF.DdsTds];
       [DispInfE1IF.DnUx,     DispInfE1IF.DssUx,    DispInfE1IF.DdsUx,    -DispInfE2IF.DnUx,    -DispInfE2IF.DssUx,    -DispInfE2IF.DdsUx   ];
       [DispInfE1IF.DnUy,     DispInfE1IF.DssUy,    DispInfE1IF.DdsUy,    -DispInfE2IF.DnUy,    -DispInfE2IF.DssUy,    -DispInfE2IF.DdsUy   ];
       [DispInfE1IF.DnUz,     DispInfE1IF.DssUz,    DispInfE1IF.DdsUz,    -DispInfE2IF.DnUz,    -DispInfE2IF.DssUz,    -DispInfE2IF.DdsUz   ];
       [-StressInfE1FB.DnTn, -StressInfE1FB.DssTn, -StressInfE1FB.DdsTn,   zerinfA,              zerinfA,               zerinfA             ];
       [-StressInfE1FB.DnTss,-StressInfE1FB.DssTss,-StressInfE1FB.DdsTss,  zerinfA,              zerinfA,               zerinfA             ];
       [-StressInfE1FB.DnTds,-StressInfE1FB.DssTds,-StressInfE1FB.DdsTds,  zerinfA,              zerinfA,               zerinfA             ];    
       [zerinfB,              zerinfB,              zerinfB,              -StressInfE2FB.DnTn,  -StressInfE2FB.DssTn,  -StressInfE2FB.DdsTn ]; %stress infs, only the influence on freeboundary
       [zerinfB,              zerinfB,              zerinfB,              -StressInfE2FB.DnTss, -StressInfE2FB.DssTss, -StressInfE2FB.DdsTss];
       [zerinfB,              zerinfB,              zerinfB,              -StressInfE2FB.DnTds, -StressInfE2FB.DssTds, -StressInfE2FB.DdsTds]];
    
    clear StressInfE1IF DispInfE1IF StressInfE1FB StressInfE2IF StressInfE2FB zerinfA zerinfB

end

 InfMatrixCheck( A );  
 %disp('inf check off in *SlipCalculator3d*, put back on'); 
    
if Option=='F'
    zer=zeros(size(DispInfE2IF.DnUz,1),1);     clear DispInfE2IF
    contvec=[zer;zer;zer;zer;zer;zer];
    B= [contvec;...
        Tn(FreeBoundary1);Tss(FreeBoundary1);Tds(FreeBoundary1);...
        Tn(FreeBoundary2);Tss(FreeBoundary2);Tds(FreeBoundary2)];
else 
    %Creating vector of tractions at every element
    B= [Tn;Tss;Tds];    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The linear equations: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Option=='C' 
    %Adding some scaling parameters (Improves Friction Solver performance). 
    %We scale by the average triangle size and the shear mod. 
    [~,HlfPerim ] = AreaOfTriangle3d( P1(:,1),P1(:,2),P1(:,3),P2(:,1),P2(:,2),P2(:,3),P3(:,1),P3(:,2),P3(:,3) );
    HlfPerim=mean(HlfPerim);
    Scl=(HlfPerim/mean(mu)); 
    A=A.*Scl; 
end


%  Solve the matrix system [A][D]=[B] to get the displacement discontinuity vector D
%If we have fixed displacements the matrix is no longer square, in this
%case \ division is faster on sparse matrices 
if any(Fdisp)==1 && isa(A,'double') %sparse doesn't work for singles
    %D = A\B; disp('sparse off, line 266 SlipCalc3d')
    D = sparse(A)\B;
    %If not we are using a dense square matrix, this is faster if square as we
    %do not need to spend time allocating this as sparse. 
else
    D = A\B;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grabbing the data back
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  Extract displacements from D vector into subvectors
if Option=='B' || Option=='E'     

    [ Dn,Dss,Dds ] = ExtractArraysFromVector( D );
    
elseif Option=='D'      
    
    [ Dss,Dds ] = ExtractArraysFromVector( D );
    Dn=zeros(size(Dss));
    
elseif Option=='F'  
    
    DE1=D(1:(3*NUME1),:);
    DE2=D((3*NUME1)+1:end,:);
    [ DnE1,DssE1,DdsE1 ] = ExtractArraysFromVector( DE1 );
    [ DnE2,DssE2,DdsE2 ] = ExtractArraysFromVector( DE2 );
    Dn  =[DnE1;DnE2];
    Dss =[DssE1;DssE2];     
    Dds =[DdsE1;DdsE2];
    
elseif Option=='C' 
    %friction solver
    [Dn,Dss,Dds] = LinearCompFrictionSolver3d(D,A,Sf,Mu,NUM,Tn,Tss,Tds);
    %Scaling back the data. 
    Dn=Dn*Scl;
    Dss=Dss*Scl;
    Dds=Dds*Scl;

end

clear A B

if Option=='E'
        Pxx=0;
        Pyy=0;
        Pzz=0;
        Pxy=0;
        Pxz=0;
        Pyz=0;        
end

end %end func

