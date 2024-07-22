% © 2022 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
classdef FirstDerivativeXFourier < FirstDerivativeOperatorFourier
    methods
        function obj = FirstDerivativeXFourier(wavenumberVector)
            arguments
                % the prescribed shape (1,:) is crucial to ensure correct broadcasting for 2D and 3D data arrays
                wavenumberVector (1,:) double {mustBeReal} = 0.0
            end
            obj = obj@FirstDerivativeOperatorFourier(wavenumberVector);
        end
        
        function derivative = compute_derivative(obj, dataArray)
            arguments
                obj FirstDerivativeXFourier
                dataArray {mustBeNumeric}  % todo: is there a way to ensure complex input?
            end
            if isvector(dataArray)  % dataArray is a vector of size 1xn or nx1: 1D data
                % make sure the leading dimension of dataArray is singleton, analogous to wavenumberVector.
                % Otherwise, the following product is interpreted as outer product, which results in a matrix
                % instead of the desired vector
                dataArray = transpose(shiftdim(dataArray));
                derivative = obj.derivativeOperator .* dataArray;
            else
                nDims = ndims(dataArray);
                if nDims == 2  % dataArray is a matrix of size mxn (with m, n not equal to 1): 2D data
                    [nz, ~] = size(dataArray);
                    derivativeOperator = repmat(obj.derivativeOperator, nz, 1);  % broadcast to 2D (DNS data layout)
                    derivative = derivativeOperator .* dataArray;
                elseif nDims == 3  % 3D data: compute x derivative at all y locations
                    [nz, ~, ny] = size(dataArray);
                    derivativeOperator = repmat(obj.derivativeOperator, nz, 1, ny); % broadcast to 3D (DNS data layout)
                    derivative = derivativeOperator .* dataArray;
                else
                    error('x derivative not implemented for size of input array.');
                end
            end
        end
    end
end