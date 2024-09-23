function [ K,E,lambda,nu,mu ] = ElasticConstantsCheck( varargin )
%ElasticConstantsCheck: Function to check the elastic constants input into
%               the function are consistent. Then allows for export of
%               constants that may be required. 
%               Input arguments must be predefined with the correct name
%               (see inputs description). This is case insensitive, the
%               inputs must also lie within reasonable elastic bounds.
% 
%
% usage #1: Get all elastic constants from two input
% [ K,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu )
%
% usage #2: Check input elastic constants match
% ElasticConstantsCheck( mu,nu,lambda )
%               
%
% Arguments: (input)
%   K           - Bulk modulus
%
%   E           - Young's modulus
%
% lambda        - Lame's constant
%
%   nu          - Poisson's ratio 
%
%   mu          - Shear modulus
%
%
% Example usage 1: Calculating other constants. 
% nu=0.25;
% mu=500; 
% [ K,E,lambda,nu,mu ] = ElasticConstantsCheck( nu,mu );
%
% Example usage 2:  Checking that the constants are compatible. 
% nu=0.25;
% mu=500; 
% lambda=6;
% ElasticConstantsCheck( lambda,nu,mu );
%
% Example usage 3: Check that the function returns same constants: 
% mu=50; 
% nu=0.01; 
% [ K,E,~,~,~ ] 		= ElasticConstantsCheck( nu,mu);
% [ K,~,lambda,~,~ ] 	= ElasticConstantsCheck(K,E);
% [ K,~,~,nu,~ ] 		= ElasticConstantsCheck(K,lambda);
% [ K,~,~,~,mu ] 		= ElasticConstantsCheck(K,nu);
% [ ~,E,lambda,~,~ ] 	= ElasticConstantsCheck(E,mu);
% [ ~,E,~,nu,~ ] 		= ElasticConstantsCheck(E,lambda);
% [ ~,E,~,~,mu ] 		= ElasticConstantsCheck(E,nu);
% [ ~,~,lambda,nu,~ ]   = ElasticConstantsCheck(E,mu);
% [ ~,~,lambda,~,mu ]   = ElasticConstantsCheck(lambda,nu);
% [ ~,~,~,nu2,mu2 ]       = ElasticConstantsCheck(lambda,mu);
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%% Part 1: getting inputs and checking these are within reasonable lims. 
if nargin<2
    error('more than one elastic constant must be supplied')
end

%Getting the input arguments.
VarNames=[];
for i=1:numel(varargin)
    %Grab current varargins input name
    VarNames{i}=inputname(i); 
end

%Comparing strings (Case insensitive). 
KExist      = strcmpi('K',      VarNames);
EExist      = strcmpi('E',      VarNames);
LambdaExist = strcmpi('lambda', VarNames);
nuExist     = strcmpi('nu',     VarNames);
muExist     = strcmpi('mu',     VarNames);

%Getting the elastic constants out of input args and checking for errors. 
if any(KExist)
    K=varargin{find(KExist)};
    Kinp=K; %Assigning to check it doesnt change
    if K<0
        error('Bulk modulus is negative')
    end
end
if any(EExist)
    E=varargin{find(EExist)};
    Einp=E; %Assigning to check it doesnt change
    if E<0
        error('Young s modulus is negative')
    end
end
if any(LambdaExist)
    lambda=varargin{find(LambdaExist)};
    Lambdainp=lambda; %Assigning to check it doesnt change
end
if any(nuExist)
    nu=varargin{find(nuExist)};
    nuinp=nu; %Assigning to check it doesnt change
    if nu<-1 || nu>0.5
        error('Poisson s ratio  out of limits')
    end
end
if any(muExist)
    mu=varargin{find(muExist)};
    muinp=mu; %Assigning to check it doesnt change
    if mu<0
        error('Shear modulus is negative')
    end
end


%% Part 2: Now calculating other elastic constants:

if any(KExist) && any(EExist)
    %1. K and E.
    lambda=((3*K)*((3*K)-E))/((9*K)-E);
    nu=((3*K)-E)/(6*K);
    mu=(3*K*E)/((9*K)-E);
end

if any(KExist) && any(LambdaExist)
    %2. K and lambda.
    E=((9*K)*(K-lambda))/((3*K)-lambda);
    nu=lambda/((3*K)-lambda);
    mu=(3/2)*(K-lambda);
end

if any(KExist) && any(nuExist)
    %3. K and nu.
    E=(3*K)*(1-(2*nu));
    lambda=(3*K*nu)/(1+nu);
    mu=((3*K)*(1-(2*nu)))/(2*(1+nu));
end

if any(KExist) && any(muExist)
    %4. K and mu
    E=(9*K*mu)/((3*K)+mu);
    lambda=((3*K)-(2*mu))/3; 
    nu=((3*K)-(2*mu))/((6*K)+(2*mu)); 
end

if any(EExist) && any(LambdaExist)
    %1. E and Lambda
    R=sqrt((E^2)+(9*(lambda^2))+(2*E*lambda));
    K=(E+(3*lambda)+R)/6;
    nu=(2*lambda)/(E+lambda+R);
    mu=(E-(3*lambda)+R)/4;
end

if any(EExist) && any(nuExist)
    %2. E and nu
    K=E/(3*(1-(2*nu)));
    lambda=(E*nu)/((1+nu)*(1-(2*nu)));
    mu=E/(2*(1+nu));
end

if any(EExist) && any(muExist)
    %3. E and mu
    K=(E*mu)/(3*((3*mu)-E));
    lambda=(mu*(E-(2*mu)))/((3*mu)-E);
    nu=(E-(2*mu))/(2*mu);
end

if any(LambdaExist) && any(nuExist)
    %1. Lambda and nu
    K=(lambda*(1+nu))/(3*nu);
    E=(lambda*(1+nu)*(1-(2*nu)))/nu;
    mu=(lambda*(1-(2*nu)))/(2*nu);
end

if any(LambdaExist) && any(muExist)
    %2. Lambda and mu
    K=((3*lambda)+(2*mu))/3;
    E=(mu*((3*lambda)+(2*mu)))/(lambda+mu);
    nu=lambda/(2*(lambda+mu));
end

if any(nuExist) && any(muExist)
    %1. nu and mu
    K=((2*mu)*(1+nu))/(3*(1-(2*nu)));
    E=(2*mu)*(1+nu);
    lambda=(2*mu*nu)/(1-(2*nu));
end


%% Part 3: Now checking all input constants match those calculated in func
if nargin>2
    if any(KExist)
        if Kinp~=K
            error('Input elastic constants are not consistent')
        end
    end
    
    if any(EExist)
        if Einp~=E
            error('Input elastic constants are not consistent')
        end
    end
    
    if any(LambdaExist)
        if Lambdainp~=lambda
            error('Input elastic constants are not consistent')
        end
    end
    
    if any(nuExist)
        if nuinp~=nu
            error('Input elastic constants are not consistent')
        end
    end
    
    if any(muExist)
        if muinp~=mu
            error('Input elastic constants are not consistent')
        end
    end
    
end

%Finishing func

end