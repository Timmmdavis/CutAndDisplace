function [ varargout ] = ExtractCols( Data,startat )
%ExtractCols grabs the x columns set as the number of output arguments. 
%startat is the column number you want to start at. its just a single value

%   Copyright 2017, Tim Davis, The University of Aberdeen


if nargin ==1
    startat=0; %start at col 1. 
else
    startat=startat-1; %we add start at on. Starting at row 1 we do not need to add to loop i. 
end

varargout= cell(1,nargout); %Preallocate cell array
for i=1:nargout
    i=i+startat;
    varargout{i}=Data(:,i);
end

% ColumnNo=size(Data,2); %Sz=1=(Rows down),2=(Cols across)
% chk=nargout==ColumnNo;
% if chk==0
%     error('number of aruguments out do not match the column numbers of the input')
% end    


end

