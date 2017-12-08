function PctOfData = Percentile(Data, Prctile)
% Percentile: Calculates the percentile of data without the stats toolbox.
%                     Not exactly the same result but very close.
%                     This assumes you are using this like MATLAB's
%                     'prctile' function but without the 3rd argument.
%                     Based on comment in:
%                     https://de.mathworks.com/matlabcentral/answers/127944
%                     -how-to-pick-the-j-th-percentile-of-a-vector
%
% usage #1:
% PctOfData = Percentile(Data, Prctile)
%
% Arguments: (input)
% Data              - The input data that is used to find the percentile. 
%
% Prctile           - Percentage defined between the interval of [0,100].
%
%
% Arguments: (output)
% PctOfData         - The percentile of the input data at the value defined
%                     in Prctile. The output 'PctOfData's rows are the
%                     percentiles asked for and the the cols the
%                     percentiles for each col of input 'Data'.
%
% Example usage:
%
% %As MATLAB's example:
% X = (1:5)'*(2:6);
% Y = Percentile(X,[25,50]);
% % And if you have the stats toolbox compare with
% Ymtlb = prctile(X,[25,50]);
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Making sure this deals with row vectors. 
if isvector(Data)
    if ~iscolumn(Data)
        Data=Data';
    end 
end

%Going to fill each row no with the differnt percentiles of the cols. 
%Going to fill each col no with the results of each col.
PctOfData=zeros(numel(Prctile),size(Data,2));

pctl = @(v,p) interp1(linspace(0.5/length(v), 1-0.5/length(v), length(v))', sort(v), p*0.01, 'spline');

%Do for every column
for jj=1:size(Data,2) 
    %Get first col
    arrlp=Data(:,jj);
    %Do for every percentile input
    for i=1:numel(Prctile)
        %Calc using func above. 
        PctOfData(i,jj) = pctl(arrlp, Prctile(i));
    end
end

end
