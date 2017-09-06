function [TensileSlipDisp,StrikeSlipDisp,DipSlipDisp] = LinearCompFrictionSolver3d(D,A,Sf,Mu,ne,TnDr,TssDr,TdsDr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Complementarity formulation
% Reformulate as a complementarity problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Copyright 2017, Tim Davis, The University of Aberdeen
%Inverted to be similar to C in Ritz script
C = inv(A);
                                                    clear A Atots B

% Construct M and Q for the equation f(Z)=M*Z+Q, where
%  f(Z) = |-Dn|  and  Z = |-tn|   
%         |Dss+|          |tss+|   
%         |Dds+|          |tds+|    
%         |Tss-|          |Dss-|   
%         |Tds-|          |Dds-|    
%Now extract variables from inverted matrix C % Rename each submatrix for easier assmbly into M. Equation 28 - Kaven 2012
% Taking shear 2 as ss and shear 3 as ds
CDnTn=C(1:ne,1:ne);					%Ann
CDssTn=C(1:ne,ne+1:2*ne);			%An2
CDdsTn=C(1:ne,2*ne+1:3*ne);  		%An3
CDnTss=C(ne+1:2*ne,1:ne);			%A2n 
CDssTss=C(ne+1:2*ne,ne+1:2*ne);		%A22
CDdsTss=C(ne+1:2*ne,2*ne+1:3*ne);	%A23
CDnTds=C(2*ne+1:3*ne,1:ne);			%A3n
CDssTds=C(2*ne+1:3*ne,ne+1:2*ne);	%A32
CDdsTds=C(2*ne+1:3*ne,2*ne+1:3*ne);	%A33
                                                    clear C B
%Making sure input frictions are vectors
Mu=zeros(ne,1)+Mu; 
Sf=zeros(ne,1)+Sf; 

% Form ne by ne array with coefficients of friction on the diagonal.
dMU = diag(Mu);
% Allocate ne by ne identity and zero matrices.
ID = eye(ne);
ZE = zeros(ne);
% Input frictional strength, Sf, Column vector exists.

%Original formulation as ritz 2012
% Construct matrix M. Modified form of Equation 28 - Kaven 2012
M =  [(CDnTn -(CDssTn*dMU)  -(CDdsTn*dMU)),  CDssTn,  CDdsTn,  ZE, ZE;
      (CDnTss-(CDssTss*dMU) -(CDdsTss*dMU)), CDssTss, CDdsTss, ID, ZE;
	  (CDnTds-(CDssTds*dMU) -(CDdsTds*dMU)), CDssTds, CDdsTds, ZE, ID;
      (2*dMU),                               -ID,      ZE,      ZE, ZE;
	  (2*dMU),                                ZE,     -ID,      ZE, ZE];

% Construct 3ne by 1 column vector Q.
Q = [(D(1:ne,:)              -(CDssTn*Sf) -(CDdsTn*Sf)); 
	 (D((ne+1):(2*ne),:)     -(CDssTss*Sf)-(CDdsTss*Sf));
	 (D(((2*ne)+1):(3*ne),:) -(CDssTds*Sf)-(CDdsTds*Sf));
	 2*Sf; 
	 2*Sf];

                                                    clear CDnTn CDssTn CDdsTn CDnTss CDssTss CDdsTss CDnTds CDssTds CDdsTds ZE ID
                                                    clear dMU DD Str
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                    
% Solve the complementarity problem using the PATH\LCP algorithm.
% Source: http://pages.cs.wisc.edu/~ferris/path/
% pathlcp.m must be in the same MATLAB directory.
tic;
disp('StartingLCP')
%%%%%%%path
%Z= pathlcp(M,Q);
%%%%%%%path

%%%%%%%LCP solve
Z = fischer_newton3d(M,Q); disp('PreAllocating sparse in fischer_newton would be faster')
close %clears the last fig 
%%%%%%%LCP solve

Timer = toc;
disp('LcpSpeed(s)');disp(Timer);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fZ = M*Z+Q;
                                                    clear M Q
                                                    
%Extracting slip                                      
TensileSlipDisp =  -fZ(1:ne,:);                            
StrikeSlipDisp = Z(3*ne+1:4*ne,1)-fZ(ne+1:2*ne,1);     
DipSlipDisp = Z(4*ne+1:5*ne,1)-fZ(2*ne+1:3*ne,1);

disp ('check Dss convention, different to normal solver when last checked')


%Correcting the slip convention as the functions of Nikkhoo M. as these use the opposite sign convention to Ole Kavens frictional solution for slip.  
TensileSlipDisp=-TensileSlipDisp;
StrikeSlipDisp=-StrikeSlipDisp;
DipSlipDisp=-DipSlipDisp;

% %Stress driving displacement on elements, after frictional resistance is overcome
% %You can put these into the regular C&S non frictiona; solver and get the same slip. 
TnNeg=Z(1:ne,1); %compressive stress (if it exists but with the wrong sign, positive should be neg, engin conv) 
Tn = TnDr+TnNeg; %Driving stress + compressive stress with wrong sign
Tss = TssDr+((Z(ne+1:2*ne,1)-(fZ(3*ne+1:4*ne,1)))-(Mu.*TnNeg)+Sf); % =Ts+(TsLL-DnRL)-Mu*TnNeg+Sf
Tds = TdsDr+((Z(2*ne+1:3*ne,1)-(fZ(4*ne+1:5*ne,1)))-(Mu.*TnNeg)+Sf); % =Ts+(TsLL-DnRL)-Mu*TnNeg+Sf

end

