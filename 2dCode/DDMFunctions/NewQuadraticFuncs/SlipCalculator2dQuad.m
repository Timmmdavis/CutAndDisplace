function [ShearDisp,TensileDisp,Pxx,Pyy,Pxy]=SlipCalculator2dQuad(x,y,xe,ye,a,Beta,...
Pxx,Pyy,Pxy,Tn,Ts,nu,E,halfspace,NUM,NormAng,strain,Fdisp,Mu,Sf,Option)
%SlipCalculator Calculates the 2 slip components on the imported fault line
%from  boundary condition defined by the user.
%
%   Influence coefficient functions from:
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
%   Mu, coefficient of friction
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
        [ Pxx,Pyy,Pxy ] = Hooke'sLaw2dStrain2Stress( Pxx,Pyy,Pxy,E,nu,(E/(2*(1+nu))));
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
    
    %Replicating for quad els
    [ Tn,Ts ] = RepeatQuadValues( Tn,Ts );
    
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
    [ DsTn,DsTs,DnTn,DnTs,Ds_Ux,Ds_Uy,Dn_Ux,Dn_Uy ] = CalculateInfluencematrices2dQuad(NUM,halfspace,x,y,xe,ye,a,Beta,nu,E,NormAng,Fdisp );
                                                                               
else %Inhomo calc
        %not hooked up for quadratic els
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

    %Repeating normals for quad collation pnts
    [ nx3,ny3 ] = RepeatQuadValues( nx,ny );
    
    %Making remote Sxx,Syy,Sxy a list the size of the number of elements
    ZeroList = zeros(numel(nx3),1);
    
    Sxx = ZeroList+Pxx;
    Syy = ZeroList+Pyy;
    Sxy = ZeroList+Pxy;
    
    %Calculates the cartesian components of traction on the elements
    [ Tx,Ty ] = TractionVectorCartesianComponents2d( Sxx,Syy,Sxy,nx3,ny3 );
    
    %Calculates the normal and shear tractions on the elements      
    [ Tn,Ts ] = CalculateNormalShearTraction2d( Tx,Ty,nx3,ny3);

end

    %for option E traction has been already defined
if Option=='E'
    %its already there, just replicate for quads
    [ Tn,Ts ] = RepeatQuadValues( Tn,Ts );

end

    clear ny nx ZeroList
    %clear  Pxx Pyy Pxy 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Creating large arrays from the influence matrices and the traction matrices.
    %performing the system of linear equations to find the slip due to the boundary conditions. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %Loop that forces any fixed elements to stick by changing the influnece
    %matrices for these elements into disp inf matrices. 
    [ Fdisp ] = RepeatQuadValues( Fdisp );
    
if Option~='F'    
    if any(Fdisp)==1
    [ DsTn,DsTs,DnTn,DnTs,Tn,Ts,NUM ]...
    = FixingDisp_InfRowSwap2d_Quad( DsTn,DsTs,DnTn,DnTs,Ds_Ux,Ds_Uy,Dn_Ux,Dn_Uy,NUM,Tn,Ts,Fdisp);
    clear Ds_Ux Ds_Uy Dn_Ux Dn_Uy
    end
else %We need to do this for both elastics
        %Inhomo calc is not hooked up for quadratic els
end    
    
if Option=='B' || Option=='E' || Option=='C'  
    Ats = [-DnTn,-DsTn]; %Convention of C&S inf coeffs is geological see fig 5.1 of C&S
    Ash = [-DnTs,-DsTs]; %Convention of C&S inf coeffs is geological see fig 5.1 of C&S
    clear DsTs DnTs DsTn DnTn                                              
    A=  [Ats;Ash];
    clear Ash Ats 
elseif Option=='D'                                                
    A = [-DsTn;-DsTs];
elseif Option=='F'         
        %Inhomo calc is not hooked up for quadratic els
end
   
    %checking the created influence matrix
    [ A ] = InfMatrixCheck( A ); 
    
    %Creating vector of tractions at every element
if Option=='F'
        %Inhomo calc is not hooked up for quadratic els

elseif Option=='B' || Option=='D' || Option=='E'  
    B = [Tn;Ts];  
elseif Option=='C'
    B = [Tn;Ts]; %[Tn;-Ts]; in linear 
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %The linear equations: 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

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
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Grabbing the data back
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %  Extract displacements from D vector into subvectors
  
if Option=='B' || Option=='E'     
    TensileDisp= D(1:NUM*3);          % Dn elements are in upper half of D
    ShearDisp = D((NUM*3)+1:2*(NUM*3));     % Ds elements are in lower half of D   
elseif Option=='D'      
    ShearDisp=D;                    % Ds elements in D 
    TensileDisp=zeros(size(D));     % No Dn by its nature
elseif Option=='F'  
        %Inhomo calc is not hooked up for quadratic els
elseif Option=='C' 
    %friction solver
    [ShearDisp,TensileDisp]=LinearCompFrictionSolver(D,A,Sf,Mu,NUM*3,Tn,Ts);
end  
    clear A B

if Option=='E'
    Pxx=0;
    Pyy=0;
    Pxy=0;
end
    
end %end entire func
