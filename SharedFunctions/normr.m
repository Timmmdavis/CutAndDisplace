function [ NMat ] = normr( Mat )
% normr: Normalises the rows of a matrix (vector norm). Like MATLAB's normr 
%               function but extra toolbox not needed
%
% usage #1: 
% [ NMat ] = normr( Mat )
%                 
% Arguments: (input)
% Mat           - The matrix that will have its rows normalised.
%
% Arguments: (input)
% NMat          - The input matrix but each row is normalised.
%
%
% Example usage:
%
% m = [1 2; 3 4];
% mn=normr(m);
% Firstrow=sqrt(mn(1)^2+mn(3)^2);
% disp(Firstrow);
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Length of each row
[~,NoCols]=size(Mat);

if (NoCols == 1)
    NMat=Mat./abs(Mat);
else
    NMat=sqrt(1./sum(Mat.^2,2))*ones(1,NoCols).*Mat;
end

end

