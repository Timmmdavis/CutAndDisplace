function [Dn,Dss,Dds] = LinearCompFrictionSolver3d(D,A,Sf,Mu,ne,TnDr,TssDr,TdsDr)
% LinearCompFrictionSolver3d: Formulates the results of the BEM calculation
%                     as a complementarity problem. 
%                     See:
%                     Ritz, E., Mutlu, O. and Pollard, D.D., 2012.
%                     Integrating complementarity into the 2D displacement
%                     discontinuity boundary element method to model faults
%                     and fractures with frictional contact properties.
%                     Computers & geosciences, 45, pp.304-312.
%                       
%                     And:
%
%                     Kaven, J.O., Hickman, S.H., Davatzes, N.C. and Mutlu,
%                     O., 2012. Linear complementarity formulation for 3D
%                     frictional sliding problems. Computational
%                     Geosciences, 16(3), pp.613-624.
%                   
%               
% usage #1:
% [ Ds,Dn ] = LinearCompFrictionSolver3d(D,A,Sf,Mu,ne,TnDr,TsDr)
%
% Arguments: (input)
% D                 - Element displacements calculated normally, i.e. these
%                     can interpenetrate. First 3rd of this col vector is
%                     normal disp and the second strike slip and 3rd is dip
%                     slip.
%
% A                 - Typical 3D DDM BEM influence matrix. 
%
% Sf                - Either list or single value of the cohesion (sliding
%                     friction) at every element.
%
% Mu                - Either list or single value of the coefficient of
%                     friction on every element.
%
% ne                - Number of elements.
%
% TnDr              - The normal traction at each element (remote stress). 
% 
% TssDr             - The strike-slip traction at each element (remote stress). 
%
% TdsDr             - The dip-slip traction at each element (remote stress). 
%
% Arguments: (output)
% Dn,Dss,Dds        - The normal, strike-slip and dip-slip displacement at 
%                     every element including frictional constraints. 
%
%
% Example usage:
%
% N/A
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Inverted to be similar to C in Ritz script
C = inv(A);   clear A 
% Another way to get [C], Slower but more intuitive.
% Z=zeros((NUM*3),1);
% for i=1:(NUM*3)
% Z=Z*0;
% Z(i)=1;
% C(:,i)=A\Z;
% end
% C=C';

%Therefore each col in [C] represents how much each 
%element must displace to cause a traction
%of one unit at element i.

% Construct a and b for the equation [y]=[a][x]+[b], where
%  [y]  = |+Dn|     and [x] = |-tn|   
%         |+Dss|              |-tss|   
%         |+Dds|              |-tds|    
%         |+Tss|              |-Dss|   
%         |+Tds|              |-Dds|   
% Where plus and minus follow the convention that positive (+) represents: 
% Dn = Opening
% Dss= Left lat movement
% Dds= Reverse slip when normals face up
% Tn = Tension
% Tss= Positive when counter clockwise from normal
% Tds= Positive when facing in positive z dir on tris normal side. 

% Creating some lengths we can use for extraction of sub matricies/vectors
% ne being the number of elements in one edge a submatrix/vector
L1=1:ne;
L2=ne+1:2*ne;
L3=2*ne+1:3*ne;
L4=3*ne+1:4*ne;
L5=4*ne+1:5*ne;


% Now extract variables from inverted matrix [C] 
%Submatrix in 'column' 1, disps cause traction Tn
CDnTn   =C(L1,L1);
CDssTn  =C(L2,L1);	
CDdsTn  =C(L3,L1);
%Submatrix in 'column' 2, disps cause traction Tss
CDnTss  =C(L1,L2);	
CDssTss =C(L2,L2);	
CDdsTss =C(L3,L2);
%Submatrix in 'column' 3, disps cause traction Tds
CDnTds  =C(L1,L3);  		
CDssTds =C(L2,L3);	
CDdsTds =C(L3,L3);	
                                                    clear C 
%Making sure input frictions are vectors
Mu=zeros(ne,1)+Mu; 
Sf=zeros(ne,1)+Sf; 

% Form ne by ne array with coefficients of friction on the diagonal.
dMU = diag(Mu);
% Allocate ne by ne identity and zero matrices.
ID = eye(ne);
ZE = zeros(ne);

% Construct matrix [a]. Modified form of Equation 28 - Kaven 2012
a =	[(CDnTn -(CDnTss*dMU)  -(CDnTds*dMU)),  CDnTss,  CDnTds,  ZE, ZE;
	 (CDssTn-(CDssTss*dMU) -(CDssTds*dMU)), CDssTss, CDssTds, ID, ZE;
	 (CDdsTn-(CDdsTss*dMU) -(CDdsTds*dMU)), CDdsTss, CDdsTds, ZE, ID;
	 (2*dMU),                               -ID,      ZE,     ZE, ZE;
	 (2*dMU),                                ZE,     -ID,     ZE, ZE];

% Construct column vector [b].
b =	[(D(L1,:)-(CDnTss*Sf) -(CDnTds*Sf)); 
	 (D(L2,:)-(CDssTss*Sf)-(CDssTds*Sf));
	 (D(L3,:)-(CDdsTss*Sf)-(CDdsTds*Sf));
	 2*Sf; 
	 2*Sf];
                                                    clear CDnTn CDssTn CDdsTn CDnTss CDssTss CDdsTss CDnTds CDssTds CDdsTds 
                                                    clear dMU ZE ID
                                                    clear CDnTn CDssTn CDdsTn CDnTss CDssTss CDdsTss CDnTds CDssTds CDdsTds 
                                                    clear dMU ZE ID
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                    
% Solve the linear complementarity problem (calling function that does this)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                    
tic; %Timing run time
disp('StartingLCP')

%%%%%%%path
%x= pathlcp(a,b);
%%%%%%%path

%%%%%%%LCP solve
x = fischer_newton(a,b); 
close %clears the last fig 
%%%%%%%LCP solve

Timer = toc;
disp('LcpSpeed(s)');disp(Timer);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculating [y], comp eq: [y]=[a][x]+[b]
y = a*x+b;
                                                    clear a b 
                                                    
%Extracting sub parts of vectors: 
%[x]
x1=x(L1,1);
x2=x(L2,1);
x3=x(L3,1);
x4=x(L4,1);
x5=x(L5,1);
%[y]
y1=y(L1,1);
y2=y(L2,1);
y3=y(L3,1);
y4=y(L4,1);
y5=y(L5,1);
                                                    clear L1 L2 L3 L4 L5

%Extracting slip                                      
Dn  = y1;                            
Dss = y2-x4;     
Dds = y3-x5;

%Extracting the resultant traction at each elements midpoint
Tn=-x1;
Tss=y4-x2;
Tds=y5-x3;
%This the remote stress left over once stress induced by the slip
%from the boundary has been removed

end

