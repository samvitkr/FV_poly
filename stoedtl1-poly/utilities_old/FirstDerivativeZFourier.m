% © 2022 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
classdef FirstDerivativeZFourier < FirstDerivativeOperatorFourier
    methods
        function obj = FirstDerivativeZFourier(wavenumberVector)
            arguments
                % the prescribed shape (:,1) is crucial to ensure correct broadcasting for 2D and 3D data arrays
                wavenumberVector (:,1) double {mustBeReal} = 0.0
            end
            obj = obj@FirstDerivativeOperatorFourier(wavenumberVector);
        end
        
        function derivative = compute_derivative(obj, dataArray)
            arguments
                obj FirstDerivativeZFourier
                dataArray {mustBeNumeric}  % is there a way to ensure complex input?
            end
            if isvector(dataArray)  % dataArray is a vector of size 1xn or nx1: 1D data
                % make sure the leading dimension of dataArray is non-singleton, analogous to the wavenumberVector.
                % Otherwise, the following multiplication is interpreted as outer product, which results in a matrix
                % instead of the desired vector.
                dataArray = shiftdim(dataArray);
                derivative = obj.derivativeOperator .* dataArray;
            else
                nDims = ndims(dataArray);
                if nDims == 2  % dataArray is a matrix of size mxn (with m, n not equal to 1): 2D data
                    [~, nx] = size(dataArray);
                    derivativeOperator = repmat(obj.derivativeOperator, 1, nx);  % broadcast to 2D (DNS data layout)
                    derivative = derivativeOperator .* dataArray;
                elseif nDims == 3  % 3D data: compute z derivative at all y locations
                    [~, nx, ny] = size(dataArray);
                    derivativeOperator = repmat(obj.derivativeOperator, 1, nx, ny); % broadcast to 3D (DNS data layout)
                    derivative = derivativeOperator .* dataArray;
                else
                    error('z derivative not implemented for size of input array.');
                end
            end
        end
    end
end