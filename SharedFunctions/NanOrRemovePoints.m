function [varargout] = NanOrRemovePoints( flag,varargin )
% NanOrRemovePoints: NaN's (or removes) data specified by a flag that
%                   matches the size of the varargin.
%               
% usage #1:
% [varargout] = NanOrRemovePoints( flag,varargin )
%
% Arguments: (input)
% flag            - Flag that signals data you would like to remove/Nan. 
%
% varargin        - Arrays that match the size and shape of the flag that
%                   will be operated on. 
%
% Arguments: (additional input)
%
% Type            - What you want to do to data that is marked by the flag.
%                   Type=0 values put to NAN
%                   Type=1 values removed from array (will turn into col vector)
%                   Call in inputs like: "'Type',1". This defaults to 0 if
%                   not input
%
% Arguments: (output)
% varargout       - The input varargin after being operated on by the flag.  
%
% Example usage (1):
%
% [X,Y]=meshgrid(-2:0.1:2,-2:0.1:2)
% flag=X>0;
% %NaN's X bigger than 0 in XY arrays. 
% [varargout] = NanOrRemovePoints( flag,X,Y )
%
% Example usage (2):
%
% [X,Y]=meshgrid(-2:0.1:2,-2:0.1:2);
% flag=X>0;
% %NaN's X bigger than 0 in XY arrays. 
% [X,Y] = NanOrRemovePoints( flag,X,Y,'Type',1 );
% disp(max(X))
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Putting Type as 0 if not defined
[ varargin,Type ]     = AdditionalArgsInVaragin( 'Type', varargin, 0 );

if Type~=1 && Type~=0
    error('Type input should be 0 or 1')
end
    
%Making sure this is a logical
flag=logical(flag);

for i=1:numel(varargin)
    
    Orig=varargin{i};
  
    if Type==0 %Change to NaN
        Orig(flag)=nan;
    
    
    else %Remove data         
        Orig(flag)=[];        
    end
    
    varargin{i}=Orig;
    
end


varargout=varargin;


end

