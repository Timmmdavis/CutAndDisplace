function [ ShearDisp,TensileDisp ] = LinearCompFrictionSolver(D,A,Sf,Mu,ne,TnDr,TsDr)
%Some of the text here is from "Ritz & Pollard 2012". 

%   Copyright 2017, Tim Davis, The University of Aberdeen

% Given the stress components, S, and the influence coefficients, A, 
%   we have a system of 2*ne simultaneous linear equations with
%   2*ne unknowns, which are the displacement discontinuities, D.
%       	  |D| =  |C||S|, where C = inv(A)
% Invert the array of influence coefficients, A, to find the array C.
 C = inv(A); 
                                                clear A    
% Complementarity formulation
% Reformulate as a complementarity problem.
% Construct M and Q for the equation f(Z)=M*Z+Q, where
%  f(Z) = |-Dn|  and  Z = |-tn|   --> D+(R),D-(L) are the magnitudes of slip in
%         | Ds+|          | ts+|   each direction, and S+,S- are the slack
%         | ts-|          | Ds-|   variables for each direction.
% Rename each submatrix for easier assmbly into M.
CDnTn= C(1:ne,1:ne);
CDsTn= C(1:ne,ne+1:2*ne);
CDnTs= C(ne+1:2*ne,1:ne); 
CDsTs= C(ne+1:2*ne,ne+1:2*ne);
                                                clear C2
% Form ne by ne array with coefficients of friction on the diagonal.
dMU = diag(Mu);

% Allocate ne by ne identity and zero matrices.
ID = eye(ne);
ZE = zeros(ne);

% %Original Ritz formulation:
% % Construct 3ne by 3ne matrix M. 
M = [(CDnTn-CDsTn*dMU),  CDsTn,  ZE;     
     (CDnTs-CDsTs*dMU),  CDsTs,  ID;
     (2*dMU),            -ID,    ZE];
% %Construct 3ne by 1 column vector Q.  
Q =[(D(1:ne))-(CDsTn*Sf);  
    (D(ne+1:2*ne)-(CDsTs*Sf));
    (2*Sf)]; 


                                                    clear CNN CNS dMU CSS
                                                    clear D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                    
% Solve the complementarity problem using the PATH algorithm.
% Source: http://pages.cs.wisc.edu/~ferris/path/
% pathlcp.m must be in the same MATLAB directory.
tic;
disp('StartingLCP')
%%%%%%path
%Z= pathlcp(M,Q);
%%%%%%path

% %%%%%%%LCP solve
 Z = fischer_newton2d(M,Q); 
 close %clears the last fig 
% %%%%%%%LCP solve

Timer = toc;
disp('LcpSpeed(s)');disp(Timer);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculate Ss, Sn, Ds, Dn on each element.
% Calculate complementary equation.
fZ = M*Z+Q;
                                                    clear M Q

% Re-calculate vectors Ss, Sn, Ds, Dn from fZ and Z.
%Slip on elements 
TensileDisp = -fZ(1:ne,1); % =-(-Dn)                                 
ShearDisp = Z(2*ne+1:3*ne,1)-fZ(ne+1:2*ne,1); % =DL-DR      

TensileDisp=-TensileDisp;
ShearDisp=-ShearDisp;

% %Stress driving displacement on elements, after frictional resistance is overcome
% %You can put these into the regular C&S non frictiona; solver and get the same slip. 
TnNeg=Z(1:ne,1); %compressive stress (if it exists but with the wrong sign, positive should be neg, engin conv) 
Tn = TnDr+TnNeg; %Driving stress + compressive stress with wrong sign
Ts = TsDr+((fZ(2*ne+1:3*ne,1))-(Z(ne+1:2*ne,1))-Mu.*TnNeg+-Sf); % =Ts+(TsLL-TSRL)-Mu*TnNeg+Sf

end

