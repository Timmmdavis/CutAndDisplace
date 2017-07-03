function [ InfluenceMatrix ] = InfMatrixCheck( InfluenceMatrix )
%Checks the influence matrix for nans and condition number. 
%Both introduce unacceptable errors. 

%   Copyright 2017, Tim Davis, The University of Aberdeen

bad = isnan(InfluenceMatrix);
nanbad=max(bad);
if nanbad==1
    disp 'Nan values exist in your inf matrix, you should find out how these got there. Overlapping segments?')
InfluenceMatrix(bad) = 0;
end

c = cond(InfluenceMatrix,2);
if c>5000 %I have chosen 5000 arbitarily, it just worked ok for some of my problems. 
disp 'High Matrix Condition Number, this may introduce rigid body movememnts (See Crouch and Starfield 1984 p96), To check see if displacements of points close to triangle faces with calculated slip move a similar amount. If not then you should be worried that the solution is poorly handled far from your fault surface.'
disp 'Matrix Condition Number, http://math.stackexchange.com/questions/261295/to-invert-a-matrix-condition-number-should-be-less-than-what'
%error('Bad matrix')

%Another way of checking is shown below, if the indentity matrix is 1's and
%0s this is fine. If not its poorly conditioned. 
%https://www.mathsisfun.com/algebra/matrix-inverse.html
% InfluenceMatrixInv=inv(InfluenceMatrix);
% Identity=InfluenceMatrix*InfluenceMatrixInv;

end    



end

