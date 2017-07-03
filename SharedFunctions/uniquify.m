function [vals, okFlag] = uniquify(vals)
%uniquify Ensures a given set of values is unique (adds a very small number to non-unique values)
%   [vals, okFlag] = uniquify(vals)
%
%   Some MATLAB functions (e.g., interpolation) require their inputs to be
%   a unique set of values. UNIQUIFY ensures that a given data vector is
%   unique by adding very small numbers to duplicate entries in the input.
%
%   UNIQUIFY does NOT sort the values, contrary to MATLAB's unique()
%   function - it just slightly modifies any duplicate values. The data
%   remains in its original order.
%
%   If UNIQUIFY was unable to find a unique set, the returned second arg
%   (okFlag) is false (otherwise: true).
%
%   Inputs:
%     vals - vector of values (possibly non-sorted and non-unique)
%
%     Note: UNIQUIFY normally accepts numeric (single/double) vectors.
%     ^^^^  Higher-dimentional inputs (e.g., a 2D matrix) are accepted,
%           but only the first data culumn is uniquified.
%
%   Outputs:
%     vals - vector of unique values (in original order)
%     okFlag - flag indicating whether output vals are unique or not
%
%   Examples:
%     vals = uniquify(vals);
%     vals = uniquify([1,1,1]) -1;          => [0, 2.22e-16, 4.44e-16]
%     vals = uniquify(single([1,1,1])) -1;  => [0, 1.19e-7,  2.38e-7]
%
%   Algorithm:
%     Loop over all duplicate entries, adding eps*2^iter in each iteration
%     until the entire data vector is unique.
%
%     Note: the algorithm is exponential to cap its run time at O(n*lg(N))
%     ^^^^  therefore, the deltas jump from 2^-N (eps) to intmax/2=2^30
%           (N=-23 for single, -52 for double)
%           However, in most non-pathological cases the loop runs only once
%
%   Class support:
%     single, double
%
%   Change log:
%     2007-Aug-28: Improved performance & accuracy; fixed runaway exponential deltas; fixed single/double support; added comments
%     2007-Aug-25: First version posted on <a href="http://www.mathworks.com/MATLABcentral/fileexchange/loadAuthor.do?objectType=author&mfx=1&objectId=1096533#">MathWorks File Exchange</a>
%
%   See also:
%     unique, interp1, consolidate (on the File Exchange)

% Programmed by Yair M. Altman: altmany(at)gmail.com
% $Revision: 1.4 $  $Date: 2007/08/28 02:04:32 $

  %try
      % Ensure floating-point inputs
      okFlag = isempty(vals);
      if okFlag
          return;
      else
          try
              classOk = isfloat(vals);
          catch
              classOk = isnumeric(vals);  % for MATLAB 6 compatibility
          end
          if ~classOk
              error('Input data class must be ''single'' or ''double''.');
          end
      end

      % Sort the data values to find duplicate values
      [sortedVals, sortedIdx] = sort(vals);
      diffVals = diff(sortedVals);
      badIdx = find(diffVals == 0);  % Find indices of duplicates

      % Find the minimal required delta for each of the input values
      % Note: automatically handles singles/doubles and crashes on other classes
      try
          epses = eps(vals);
      catch
          epses = repmat(eps,size(vals));  % for MATLAB 6 compatibility
      end

      % The number of loop iterations depends on the minimal step-size
      maxIters = 30 - log2(min(epses));  % intmax = 2^31-1

      % Main processing loop (normally runs just once)
      iter = 0;
      while ~isempty(badIdx) & (iter < maxIters)  %#ok for MATLAB 6 compatibility

          % Compute and add delta(s) to next higher value(s)
          deltas = 2^iter * epses(sortedIdx(badIdx)) .* badIdx;
          vals(sortedIdx(badIdx+1)) = vals(sortedIdx(badIdx+1)) + deltas;

          % Recheck modified values for duplicates (should not normally be any)
          diffVals = diff(vals(sortedIdx));
          badIdx = find(diffVals == 0);
          iter = iter + 1;
      end

      % Indicate whether the unification was successful or not
      okFlag = isempty(badIdx);

  %catch
  %    handleError;
  %end

  return;  % debug point
