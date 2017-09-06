function [ShearDisp,TensileDisp,Pxx,Pyy,Pxy]=SlipCalculator2d(x,y,xe,ye,a,Beta,...
Pxx,Pyy,Pxy,Tn,Ts,nu,E,halfspace,NUM,NormAng,strain,Fdisp,Mu,Sf,Option)
%SlipCalculator Calculates the 2 slip components on the imported fault line
%from  boundary condition defined by the user.
%
%   Influence Coefficient functions from:
%   "Boundary element methods in Solid Mechanics"
%   by S.L. Crouch and A.M. Starfield, 1983
%   Big thanks to Steve Martel for putting example MATLAB scripts online
%
%   INPUTS
%   x & y are element midpoints
%   xe & ye are also element midpoints, the coefficients function looks at
%   how each xe affects each x etc. 
%   a is an array of each each elements half length
%   Beta is the orientation of each element relative to the X axis in the
%   C&S convention (See C&S section 2.8 p22 for Beta definition)
%   Pxx, Pyy & Pxy are the stress inputs
%   nu is the Poisson's ratio
%   E is the Young's modulus
%   mu is the shear modulus
%   halfspace defines if we work out the coefficientsin a half or whole
%   space
%   NUM is the number of elements
%   NormAng is now the outward normal measured from Xaxis counter
%   clockwise to the normal
%   Pxx,Pyy,Pxy,Tn,Ts are the defined stresses or tractions at the
%   elements, either defined for each element or just a single value that
%   will then be added to all elements.
%   strain is a flag that means remote stress are actually strains, these
%   are converted to stresses before this continues. 
%   Fdisp, flag of elements that have fixed displacements. Displacement
%   coefficients are created if this exists and the system of linear
%   equations changes. 
%   Mu, Coefficient of friction
%   Sf, Sliding friction
%   Option - which option the user is using. 
    %B=under remote stress, fault can interpenetrate
    %C=friction solver as in ritz 2012
    %D=remote stress, opening components not calculated. Fault satisfys
    %boundary conditions by shearing only.
    %E=tractions defined on body, Ie pressure in a magma chamber
    %F=Inhomogeneous
    
%   Copyright 2017, Tim Davis, The University of Aberdeen
%   OUTPUTS
%   ShearDisp & TensileDisp are the displacement components on each element
%   Pxx, Pyy, and Pxy are the remote boundary far field stresses 



% If Defined slip is as a uniform remote strain, checking options that have
% a uniform 'stress' defined
if Option=='B' || Option=='C' || Option=='D' || Option=='F'

     %for the imported stress making sure these are col vectors
    [ Pxx ] = RowVecToCol( Pxx );
    [ Pyy ] = RowVecToCol( Pyy );
    [ Pxy ] = RowVecToCol( Pxy );

    if strain==1
        %Hooke's law strain to stress
        [ Pxx,Pyy,Pxy ] = HookesLaw2dStrain2Stress( Pxx,Pyy,Pxy,E,nu,(E/(2*(1+nu))));
    end  
end

%Function to check the boundary stress tensors related to Z(y in 2d) components are 0 at the
%halfspace surface
if halfspace==1
    HalfSpaceBoundaryConditionsCheck2d( ye,Pyy,Pxy )
else
end

    %Checking we can overcome the frictional resistance on the fault. No
    %point waiting for this otherwise
if Option=='C'     
    
    nx=cos(NormAng); ny=sin(NormAng);    
    %Calculates the normal and shear tractions on the elements 
    [ Tx,Ty ] = TractionVectorCartesianComponents2d( Pxx,Pyy,Pxy,nx,ny );
    [ Tn,Ts ] = CalculateNormalShearTraction2d( Tx,Ty,nx,ny);
    
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
    [ DsTn,DsTs,DnTn,DnTs,Ds_Ux,Ds_Uy,Dn_Ux,Dn_Uy ] = CalculateInfluenceMatrices2d(NUM,halfspace,x,y,xe,ye,a,Beta,nu,E,NormAng,Fdisp );
    
else %Inhomo calc

    [DsTnE1FB,DsTsE1FB,DnTnE1FB,DnTsE1FB,DsTnE1IF,DsTsE1IF,DnTnE1IF,DnTsE1IF,...
    Ds_UxE1FB,Ds_UyE1FB,Dn_UxE1FB,Dn_UyE1FB,Ds_UxE1IF,Ds_UyE1IF,Dn_UxE1IF,Dn_UyE1IF,...
    DsTnE2FB,DsTsE2FB,DnTnE2FB,DnTsE2FB,DsTnE2IF,DsTsE2IF,DnTnE2IF,DnTsE2IF,...
    Ds_UxE2FB,Ds_UyE2FB,Dn_UxE2FB,Dn_UyE2FB,Ds_UxE2IF,Ds_UyE2IF,Dn_UxE2IF,Dn_UyE2IF,...
    NUME1,NUME2,nuE1,EE1,nuE2,EE2,...
    FreeBoundary1,FreeBoundary2,FdispE1,FdispE2,FreeBoundaries,FixedDisps]...
    = InhomogeneousInfluenceBuilder2d...
    (x,y,xe,ye,a,Beta,nu,E,halfspace,NormAng,Fdisp);
      
end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Traction influence matrices should now exist. 
    %Now the defined boundary conditions are converted to traction at each elements midpoint. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Option=='B' || Option=='C' || Option=='D' || Option=='F' 
    % Computing the two direction cosines of the normal of each element of the
    % fault. sin(nx) could be used instead of ny but for consistency with 3d
    % code ny is used. 
    
    if Option=='F'     
        ax=NormAng(FreeBoundaries);  
    else    
        ax=NormAng;
    end    
    nx=cos(ax);
    ny=cos((pi/2)-ax);

    %Making remote Sxx,Syy,Sxy a list the size of the number of elements
    ZeroList = zeros(numel(nx),1);
    
    Sxx = ZeroList+Pxx;
    Syy = ZeroList+Pyy;
    Sxy = ZeroList+Pxy;
    
    %Calculates the cartesian components of traction on the elements
    [ Tx,Ty ] = TractionVectorCartesianComponents2d( Sxx,Syy,Sxy,nx,ny );
    
    %Calculates the normal and shear tractions on the elements      
    [ Tn,Ts ] = CalculateNormalShearTraction2d( Tx,Ty,nx,ny);

end

    %for option E traction has been already defined
if Option=='E'
    %do nothing
end

    clear ny nx ZeroList
    %clear  Pxx Pyy Pxy 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Creating large arrays from the influence matrices and the traction matrices.
    %performing the system of linear equations to find the slip due to the boundary conditions. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %Loop that forces any fixed elements to stick by changing the influnece
    %matrices for these elements into disp inf matrices. 
if Option~='F'    
    if any(Fdisp)==1
    [ DsTn,DsTs,DnTn,DnTs,Tn,Ts,NUM ]...
    = FixingDisp_InfRowSwap2d( DsTn,DsTs,DnTn,DnTs,Ds_Ux,Ds_Uy,Dn_Ux,Dn_Uy,NUM,Tn,Ts,Fdisp);
    clear Ds_Ux Ds_Uy Dn_Ux Dn_Uy
    end
else %We need to do this for both elastics
    %E1  %%%%%%%%
    if any(FdispE1)==1
    [ DsTnE1FB,DsTsE1FB,DnTnE1FB,DnTsE1FB,null,null,NUME1 ]...
    = FixingDisp_InfRowSwap2d( DsTnE1FB,DsTsE1FB,DnTnE1FB,DnTsE1FB,Ds_UxE1FB,Ds_UyE1FB,Dn_UxE1FB,Dn_UyE1FB,NUME1,[],[],FdispE1);
    % % %Now removing any fixed disp element cols from the interface elements
    [DsTnE1IF,DsTsE1IF,DnTnE1IF,DnTsE1IF] = ExtractData(FdispE1,4,DsTnE1IF,DsTsE1IF,DnTnE1IF,DnTsE1IF );
    [Ds_UxE1IF,Ds_UyE1IF,Dn_UxE1IF,Dn_UyE1IF] = ExtractData(FdispE1,4,Ds_UxE1IF,Ds_UyE1IF,Dn_UxE1IF,Dn_UyE1IF);  
    end
    
    
    %E2 %%%%%%%% 
    if any(FdispE2)==1
    %Flipping matrices so the fixed displacement parts are at the top of the list  
    [DsTnE2FB,DsTsE2FB,DnTnE2FB,DnTsE2FB,Ds_UxE2FB,Ds_UyE2FB,Dn_UxE2FB,Dn_UyE2FB,FdispE2]...
        =FlipColVecs(DsTnE2FB,DsTsE2FB,DnTnE2FB,DnTsE2FB,Ds_UxE2FB,Ds_UyE2FB,Dn_UxE2FB,Dn_UyE2FB,FdispE2);
    [ DsTnE2FB,DsTsE2FB,DnTnE2FB,DnTsE2FB,null,null,NUME2 ]...
    = FixingDisp_InfRowSwap2d( DsTnE2FB,DsTsE2FB,DnTnE2FB,DnTsE2FB,Ds_UxE2FB,Ds_UyE2FB,Dn_UxE2FB,Dn_UyE2FB,NUME2,[],[],FdispE2); 
    %Flipping these matrices back  
    [DsTnE2FB,DsTsE2FB,DnTnE2FB,DnTsE2FB,Ds_UxE2FB,Ds_UyE2FB,Dn_UxE2FB,Dn_UyE2FB,FdispE2]...
        =FlipColVecs(DsTnE2FB,DsTsE2FB,DnTnE2FB,DnTsE2FB,Ds_UxE2FB,Ds_UyE2FB,Dn_UxE2FB,Dn_UyE2FB,FdispE2);
    % % %Now removing any fixed disp element cols from the interface elements
    [DsTnE2IF,DsTsE2IF,DnTnE2IF,DnTsE2IF] = ExtractData(FdispE2,4,DsTnE2IF,DsTsE2IF,DnTnE2IF,DnTsE2IF);
    [Ds_UxE2IF,Ds_UyE2IF,Dn_UxE2IF,Dn_UyE2IF] = ExtractData(FdispE2,4,Ds_UxE2IF,Ds_UyE2IF,Dn_UxE2IF,Dn_UyE2IF); 
    end   
    %Need to do this for Tn and Ts 
    FD=logical(FixedDisps);    
    Tn(FD)=0;
    Ts(FD)=0;
end    
    
if Option=='B' || Option=='E'     
    Ats = [-DnTn,-DsTn];
    Ash = [-DnTs,-DsTs];
    clear DsTs DnTs DsTn DnTn                                              
    A=  [Ats;Ash];
    clear Ash Ats 
elseif Option=='C'      
    Ats = [-DnTn,-DsTn]; 
    Ash = [-DnTs,-DsTs]; 
    clear DsTs DnTs DsTn DnTn                                              
    A=  [Ats;Ash];
    clear Ash Ats 
elseif Option=='D'                                                
    A = [-DsTn;-DsTs];
elseif Option=='F'         
    zerinfA=zeros(size(DsTnE1FB,1),size(Dn_UyE2IF,2));%Sz=1=(Rows down),2=(Cols across)
    zerinfB=zeros(size(DnTnE2FB,1),size(DsTsE1FB,2)); 
    %Setting up the matrix
    A=   [[DnTnE1IF,    DsTnE1IF,   -DnTnE2IF,    -DsTnE2IF]; %stress infs that must be satisfied at interface
          [DnTsE1IF,    DsTsE1IF,   -DnTsE2IF,    -DsTsE2IF];
          [Dn_UxE1IF,   Ds_UxE1IF,  -Dn_UxE2IF,  -Ds_UxE2IF]; %disp infs that must be satisfied at interface
          [Dn_UyE1IF,   Ds_UyE1IF,  -Dn_UyE2IF,  -Ds_UyE2IF];
          [DnTnE1FB,    DsTnE1FB,      zerinfA,     zerinfA]; %stress infs, only the influence on freeboundary
          [DnTsE1FB,    DsTsE1FB,      zerinfA,     zerinfA]
          [zerinfB,     zerinfB,     DnTnE2FB,    DsTnE2FB]; %stress infs, only the influence on freeboundary
          [zerinfB,     zerinfB,     DnTsE2FB,    DsTsE2FB]];      
    clear DsTnE1IF DsTsE1IF DnTnE1IF DnTsE1IF DsTnE2IF DsTsE2IF DnTnE2IF DnTsE2IF zerinfA
    clear Ds_UxE1IFDs_UyE1IFDn_UxE1IFDn_UyE1IF Ds_UxE2IFDs_UyE2IFDn_UxE2IFDn_UyE2IF
    clear DsTnE1FB DsTsE1FB DnTnE1FB DnTsE1FB
    clear zerinfB DnTnE2FB DnTsE2FB DsTnE2FB DsTsE2FB zerinfB E2FB    
end
   
    %checking the created influence matrix
    [ A ] = InfMatrixCheck( A ); 
    
    %Creating vector of tractions at every element
if Option=='F'
    zer=zeros(size(Dn_UyE2IF,1),1);
    contvec=[zer;zer;zer;zer];
    B=-[contvec;Tn(FreeBoundary1);Ts(FreeBoundary1);Tn(FreeBoundary2);Ts(FreeBoundary2)];

elseif Option=='B' || Option=='D' || Option=='E'  
    B = [Tn;Ts];  
elseif Option=='C'
    B = [Tn;Ts]; 
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %The linear equations: 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tic

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
    
    toc

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Grabbing the data back
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %  Extract displacements from D vector into subvectors
  
if Option=='B' || Option=='E'     
    TensileDisp= D(1:NUM);          % Dn elements are in upper half of D
    ShearDisp = D(NUM+1:2*NUM);     % Ds elements are in lower half of D   
elseif Option=='D'      
    ShearDisp=D;                    % Ds elements in D 
    TensileDisp=zeros(size(D));     % No Dn by its nature
elseif Option=='F'  
    TensileDispE1       = D(1:NUME1,:);                            %Elastic1
    ShearDispE1         = D(NUME1+1:2*NUME1,:);                    %Elastic1
    TensileDispE2       = D(2*NUME1+1:2*NUME1+NUME2,:);            %Elastic2
    ShearDispE2         = D((((2*NUME1)+NUME2+1)):end,:);          %Elastic2
    TensileDisp=[TensileDispE1;TensileDispE2];
    ShearDisp=[ShearDispE1;ShearDispE2];         
elseif Option=='C' 
    %friction solver
    [ShearDisp,TensileDisp]=LinearCompFrictionSolver(D,A,Sf,Mu,NUM,Tn,Ts);
end  
    clear A B

if Option=='E'
    Pxx=0;
    Pyy=0;
    Pxy=0;
end
end %end entire func
