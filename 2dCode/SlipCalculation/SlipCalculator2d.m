function [Ds,Dn,Pxx,Pyy,Pxy]=SlipCalculator2d...
    (MidPoint,HalfLength,Pxx,Pyy,Pxy,Tn,Ts,nu,E,halfspace,LineNormalVector,strain,Fdisp,Mu,Sf,Option)
% SlipCalculator2d: Calculates the 2 slip components on the
%               imported fault line from  boundary condition defined by the
%               user.
%
% usage #1:
% [Ds,Dn,Pxx,Pyy,Pxy]=SlipCalculator2d...
%     (MidPoint,HalfLength,Pxx,Pyy,Pxy,Tn,Ts,nu,E,halfspace,LineNormalVector,strain,Fdisp,Mu,Sf,Option)
%
% Arguments: (input)
%  MidPoint   - The element midpoints in X and Y. 
%
%       a     - An array of each each elements half length.
%
% Pxx,Pyy,Pxy - Remote stresses defined by user.
%
%     Tn,Ts   - Tractions on the surface defined by the user. 
%
%       nu    - The Poisson's ratio.
%
%       E     - The Young's modulus.
%
%  halfspace  - Defines if we work out the coefficientsin a half or whole
%              space.
%
% LineNormalVector - The direction cosines, CosAx (Nx) and CosAy in a list
%                   for each element. 
%
%    Strain   - Flag to say if the input stresses are actually a strain
%              tensor. If so we convert these to stresses before
%              continuing.
%
% Fdisp       - Flag telling the user if any elements are going to be fixed
%              (if this is the case displacement influence matricies are
%              required). In the case of inhomogeneous materials this is
%              more complex, see function
%              'InhomogeneousInfluenceBuilder2d.m'.
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
%       Ds,Dn - Vectors that describe how much the elements displace in the
%               normal (Dn) and shear (Ds) directions on the elements. 
%
% Pxx,Pyy,Pxy - Remote stress if it existed on import (useful if the
%              user defined strain boundary conditions. 
%
% Example usage:
% %Assuming the elastic constants and line are defined:
%  halfspace=0;
%  [Sxx,Syy,Sxy,Tn,Ts,Mu,Sf,strain ] = CreateBlankVars;
%  strain=0; 
%  Sxx = 0; 					
%  Syy = 0;
%  Sxy = 1;                  
%  Option='B'; 
% [Ds,Dn,Pxx,Pyy,Pxy]=SlipCalculator2d...
%     (MidPoint,HalfLength,Pxx,Pyy,Pxy,Tn,Ts,nu,E,halfspace,LineNormalVector,strain,Fdisp,Mu,Sf,Option)
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

NUM=numel(MidPoint(:,1));

% If Defined slip is as a uniform remote strain, checking options that have
% a uniform 'stress' defined
if Option=='B' || Option=='C' || Option=='D' || Option=='F'

     %for the imported stress making sure these are col vectors
    [ Pxx,Pyy,Pxy ] = RowVecToCol( Pxx,Pyy,Pxy );

    if strain==1
        %Hooke's law strain to stress
        [ Pxx,Pyy,Pxy ] = HookesLaw2dStrain2Stress( Pxx,Pyy,Pxy,E,nu,(E/(2*(1+nu))));
    end  
end

%Function to check the boundary stress tensors related to Z(y in 2d) components are 0 at the
%halfspace surface
if halfspace==1
    HalfSpaceBoundaryConditionsCheck(MidPoint(:,2),Pyy,Pxy)
end

%Checking we can overcome the frictional resistance on the fault. No
%point waiting for this otherwise
if Option=='C'     
    
    CosAx=LineNormalVector(:,1);
    CosAy=LineNormalVector(:,2);  
    %Calculates the normal and shear tractions on the elements 
    [ Tx,Ty ] = TractionVectorCartesianComponents2d( Pxx,Pyy,Pxy,CosAx,CosAy );
    [ Tn,Ts ] = CalculateNormalShearTraction2d( Tx,Ty,CosAx,CosAy);
    
    FricRes=Tn.*Mu; %Negative is compression
    if all(-abs(Ts)>FricRes) 
        disp('Your frictional fault will not slip, quit? (ctrl+c)')
        pause        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating influence matrices using a function. Big square matrices of
%how much every element slipping by 1 effects traction on every other element
%If displacement if fixed on any elements we also create disp influence matrices 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if Option~='F' %not inhomo calc
    
    [StressInf,DispInf]=CalculateInfluenceMatrices2d(halfspace,MidPoint,HalfLength,nu,E,LineNormalVector,Fdisp );
                                                        
else %Inhomo calc

    [StressInfE1FB,StressInfE1IF,DispInfE1FB,DispInfE1IF,StressInfE2FB,StressInfE2IF,DispInfE2FB,DispInfE2IF,...
     NUME1,NUME2,FreeBoundary1,FreeBoundary2,FdispE1,FdispE2,FreeBoundaries,FixedDisps]...
     = InhomogeneousInfluenceBuilder2d...
     (MidPoint,HalfLength,nu,E,halfspace,LineNormalVector,Fdisp);
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Traction influence matrices should now exist. 
%Now the defined boundary conditions are converted to traction at each elements midpoint. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Option=='B' || Option=='C' || Option=='D' || Option=='F' 
    
    if Option=='F'     
        %Direction cosines
        CosAx=LineNormalVector(FreeBoundaries,1);
        CosAy=LineNormalVector(FreeBoundaries,2);  
    else    
        %Direction cosines
        CosAx=LineNormalVector(:,1);
        CosAy=LineNormalVector(:,2);  
    end    

    %Making remote Sxx,Syy,Sxy a list the size of the number of elements
    ZeroList = zeros(numel(CosAx),1);
    
    Sxx = ZeroList+Pxx;
    Syy = ZeroList+Pyy;
    Sxy = ZeroList+Pxy;
    
    %Calculates the cartesian components of traction on the elements
    [ Tx,Ty ] = TractionVectorCartesianComponents2d( Sxx,Syy,Sxy,CosAx,CosAy );
    
    %Calculates the normal and shear tractions on the elements      
    [ Tn,Ts ] = CalculateNormalShearTraction2d( Tx,Ty,CosAx,CosAy);

end

    %for option E traction has been already defined
if Option=='E'
    %do nothing
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating large arrays from the influence matrices and the traction matrices.
%performing the system of linear equations to find the slip due to the boundary conditions. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Loop that forces any fixed elements to stick by changing the influnece
%matrices for these elements into disp inf matrices. 
if Option~='F'    
    
    if any(Fdisp)==1
        [StressInf,Tn,Ts,NUM]=FixingDisp_InfRowSwap2d(StressInf,DispInf,NUM,Tn,Ts,Fdisp);
        clear DispInf
    end
    
else %We need to do this for both elastics
    
    %E1  %%%%%%%%
    if any(FdispE1)==1
        [StressInfE1FB,~,~,NUME1]...
        = FixingDisp_InfRowSwap2d(StressInfE1FB,DispInfE1FB,NUME1,[],[],FdispE1);
        %Now removing any fixed disp element cols from the interface elements
        [StressInfE1IF] = ExtractData(FdispE1,2,StressInfE1IF,'Columns',1 );
        [DispInfE1IF] = ExtractData(FdispE1,2,DispInfE1IF,'Columns',1);  
    end
    
    %E2 
    if any(FdispE2)==1
        %Flipping matrices so the fixed displacement parts are at the top of the list   
        [StressInfE2FB,DispInfE2FB,FdispE2]=FlipColVecs(StressInfE2FB,DispInfE2FB,FdispE2);
        [StressInfE2FB,~,~,~]= FixingDisp_InfRowSwap2d(StressInfE2FB,DispInfE2FB,NUME2,[],[],FdispE2);
        %Flipping these matrices back  
        [StressInfE2FB,DispInfE2FB,FdispE2]=FlipColVecs(StressInfE2FB,DispInfE2FB,FdispE2);
        %Now removing any fixed disp element cols from the interface elements
        [StressInfE2IF] = ExtractData(FdispE2,2,StressInfE2IF,'Columns',1 );
        [DispInfE2IF] = ExtractData(FdispE2,2,DispInfE2IF,'Columns',1);  
    end   
    FD=logical(FixedDisps);    
    Tn(FD)=0;
    Ts(FD)=0;
end    
    
if Option=='B' || Option=='E' || Option=='C'   
  
    Ats = [-StressInf.DnTn,-StressInf.DsTn];            StressInf = rmfield(StressInf,{'DnTn','DsTn'});   
    Ash = [-StressInf.DnTs,-StressInf.DsTs];            
    clear DsTs DnTs DsTn DnTn                                              
    A=  [Ats;Ash];
    clear Ash Ats 
	
elseif Option=='D'      
                                          
    A = [-StressInf.DsTn;-StressInf.DsTs];
	
elseif Option=='F'         

    zerinfA=zeros(size(StressInfE1FB.DsTn,1),size(DispInfE2IF.DnUy,2));%Sz=1=(Rows down),2=(Cols across)
    zerinfB=zeros(size(StressInfE2FB.DnTn,1),size(StressInfE1FB.DsTs,2)); 
	
    %Setting up the matrix
    A=   [[StressInfE1IF.DnTn,  StressInfE1IF.DsTn, -StressInfE2IF.DnTn, -StressInfE2IF.DsTn ]; %stress infs that must be satisfied at interface
          [StressInfE1IF.DnTs,  StressInfE1IF.DsTs, -StressInfE2IF.DnTs, -StressInfE2IF.DsTs ];
          [DispInfE1IF.DnUx,    DispInfE1IF.DsUx,   -DispInfE2IF.DnUx,   -DispInfE2IF.DsUx   ]; %disp infs that must be satisfied at interface
          [DispInfE1IF.DnUy,    DispInfE1IF.DsUy,   -DispInfE2IF.DnUy,   -DispInfE2IF.DsUy   ];
          [-StressInfE1FB.DnTn,-StressInfE1FB.DsTn,  zerinfA,             zerinfA            ]; %stress infs, only the influence on freeboundary
          [-StressInfE1FB.DnTs,-StressInfE1FB.DsTs,  zerinfA,             zerinfA            ];
          [zerinfB,             zerinfB,            -StressInfE2FB.DnTn, -StressInfE2FB.DsTn ]; %stress infs, only the influence on freeboundary
          [zerinfB,             zerinfB,            -StressInfE2FB.DnTs, -StressInfE2FB.DsTs ]];      
		  
    clear StressInfE1IF StressInfE2IF  DispInfE1IF StressInfE1FB StressInfE2FB
	
end
   
    InfMatrixCheck( A );  
     %disp('inf check off in *slipcalculator3d*, put back on'); 
    
    %Creating vector of tractions at every element
if Option=='F'
    zer=zeros(size(DispInfE2IF.DsUy,1),1); clear DispInfE2IF
    contvec=[zer;zer;zer;zer];
    B= [contvec;Tn(FreeBoundary1);Ts(FreeBoundary1);Tn(FreeBoundary2);Ts(FreeBoundary2)];
else
    B= [Tn;Ts];  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The linear equations: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Option=='C' 
    %Adding some scaling parameters (Improves Friction Solver performance). 
    %We scale by the average line size and the youngs mod. 
    HlfLng=mean(HalfLength);
    mu=mean(E)/(2*(1+mean(nu)));
    Scl=(HlfLng./mu).*100; 
    A=A.*Scl; 
end

%  Solve the matrix system [A][D]=[B] to get the displacement discontinuity vector D
%If we have fixed displacements the matrix is no longer square, in this
%case \ division is faster on sparse matrices 
if any(Fdisp)==1 && isa(A,'double')
    D = sparse(A)\B;
    %If not we are using a dense square matrix, this is faster if square as we
    %do not need to spend time allocating this as sparse. 
else
    D = A\B;
end
    
%toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grabbing the data back
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Extract displacements from D vector into subvectors
  
if Option=='B' || Option=='E'  
    
    [ Dn,Ds ] = ExtractArraysFromVector( D );

elseif Option=='D'    
    
    [ Ds ] = ExtractArraysFromVector( D );
    Dn=zeros(size(D));     % No Dn by its nature
    
elseif Option=='F'  
    
    DE1=D(1:(2*NUME1),:);
    DnE2=D((2*NUME1)+1:end,:);
    [ DnE1,DsE1 ] = ExtractArraysFromVector( DE1 );
    [ DnE2,DsE2 ] = ExtractArraysFromVector( DnE2 );
    Dn=[DnE1;DnE2];
    Ds=[DsE1;DsE2];   
    
elseif Option=='C' 
    %friction solver
    [Ds,Dn]=LinearCompFrictionSolver2d(D,A,Sf,Mu,NUM,Tn,Ts);
	%Scaling back the data. 
    Dn=Dn*Scl;
    Ds=Ds*Scl;
end  
clear A B

if Option=='E'
	Pxx=0;
    Pyy=0;
    Pxy=0;
end

end %end entire func
