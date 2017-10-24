function [ varargout ] = ExtractCols( Data )
%ExtractCols grabs the x columns set as the number of output arguments. 

%   Copyright 2017, Tim Davis, The University of Aberdeen

varargout= cell(1,nargout); %Preallocate cell array
    for i=1:nargout
        varargout{i}=Data(:,1);
        Data=Data(:,2:end); %Remove the first column (this is more memory efficient) 
    end

end

