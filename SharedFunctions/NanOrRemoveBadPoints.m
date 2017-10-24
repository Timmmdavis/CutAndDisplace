function [varargout] = NanOrRemoveBadPoints( flag,type,varargin )
%NanOrRemoveBadPoints 
%Create a flag and send this in, perform on the inputs 'varargin' 
%varargin the arguments you are passing in. 

%   Copyright 2017, Tim Davis, The University of Potsdam

%type=0 values put to NAN
%type=1 values removed from array (shapes may change)
    
for i=1:nargin-3
    
    Orig=varargin{i};
   
    if type==0 %change to nan
    Orig(flag)=nan;
    else       %remove data  
    Orig(flag)=[];        
    end
    varargin{i}=Orig;
end


varargout=varargin;


end

