function [ ShearDisp,TensileDisp ] = LinearCompFrictionSolver(D,A,Sf,Mu,ne,TnDr,TsDr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Complementarity formulation
% Reformulate as a complementarity problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Copyright 2017, Tim Davis, The University of Aberdeen
%Ritz ref here is:
%Ritz, E., Mutlu, O. and Pollard, D.D., 2012. Integrating complementarity
%into the 2D displacement discontinuity boundary element method to model
%faults and fractures with frictional contact properties. Computers &
%geosciences, 45, pp.304-312.

%Inverted to be similar to C in Ritz script
C = inv(A);
                                                    clear A 
% Another way to get [C], Slower but more intuitive.
% Z=zeros((NUM*2),1);
% for i=1:(NUM*2)
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
%         |+Ds|               |-ts|   
%         |+Ts|               |-Ds|   

% Where plus and minus follow the convention that positive (+) represents: 
% Dn = Opening
% Ds= Left lat movement
% Tn = Tension
% Ts= Positive when counter clockwise from normal

% Creating some lengths we can use for extraction of sub matricies/vectors
% ne being the number of elements in one edge a submatrix/vector
L1=1:ne;
L2=ne+1:2*ne;
L3=2*ne+1:3*ne;

% Now extract variables from inverted matrix [C] 
%Submatrix in 'column' 1, disps cause traction Tn
CDnTn= C(L1,L1);
CDsTn= C(L2,L1);
%Submatrix in 'column' 2, disps cause traction Ts
CDnTs= C(L1,L2);
CDsTs= C(L2,L2);
                                                clear C

%Making sure input frictions are vectors
Mu=zeros(ne,1)+Mu; 
Sf=zeros(ne,1)+Sf; 

% Form ne by ne array with coefficients of friction on the diagonal.
dMU = diag(Mu);
% Allocate ne by ne identity and zero matrices.
ID = eye(ne);
ZE = zeros(ne);

% % Construct matrix [a] 
a = [(CDnTn-CDnTs*dMU),  CDnTs,  ZE;     
     (CDsTn-CDsTs*dMU),  CDsTs,  ID;
     (2*dMU),            -ID,    ZE];
 
% %Construct 3ne by 1 column vector Q.  
b =[(D(L1)-(CDnTs*Sf));  
    (D(L2)-(CDsTs*Sf));
    (2*Sf)]; 


                                                    clear CDnTn CDnTs CDsTn CDsTs
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
x = fischer_newton2d(a,b); 
%close %clears the last fig 
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
%[y]
y1=y(L1,1);
y2=y(L2,1);
y3=y(L3,1);

                                                    clear L1 L2 L3

%Extracting slip                                      
TensileDisp = y1;                            
ShearDisp  = y2-x3;      

%Extracting the resultant traction at each elements midpoint
Tn=-x1;
Ts=y3-x2;
%This the remote stress left over once stress induced by the slip
%from the boundary has been removed

end

