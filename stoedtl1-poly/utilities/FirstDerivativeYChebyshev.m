% © 2022 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
classdef FirstDerivativeYChebyshev < FirstDerivativeOperator
    properties (Access = private)
        yGridChebyshev;
    end
    
    methods
        function obj = FirstDerivativeYChebyshev(nPointsY)
            arguments
                nPointsY (1,1) {mustBeInteger, mustBeNonnegative} = 0
            end
            if nPointsY > 0
                [obj.yGridChebyshev, obj.derivativeOperator] = chebdif(nPointsY, 1);
            else
                obj.yGridChebyshev = 0;
                obj.derivativeOperator = 0;
            end
        end
        
        function derivative = compute_derivative(obj, dataArray)
            arguments
                obj FirstDerivativeYChebyshev
                dataArray {mustBeNumeric}
            end
            if isvector(dataArray)  % dataArray is a vector of size 1xn or nx1: 1D data
                % make sure the leading dimension of dataArray is non-singleton, so that the matrix-vector
                % multiplication yields the correct result
                dataArray = shiftdim(dataArray);
                derivative = obj.derivativeOperator * dataArray;
            else
                if ndims(dataArray) == 3  % 3D data: compute y derivative at all x and z locations
                    idxOrder = [3, 2, 1];
                    dataArrayPermuted = permute(dataArray, idxOrder);  % reorder data: [z, x, y] -> [y, x, z]
                    derivativePermuted = pagemtimes(obj.derivativeOperator, dataArrayPermuted);  % apply y-derivative at each x, z
                    derivative = ipermute(derivativePermuted, idxOrder);  % back to original order: [y, x, z] -> [z, x, y]
                else
                    error('y derivative not implemented for size of input array.');
                end
            end
        end
        
        function yGrid = get_chebyshev_grid(obj)
            yGrid = obj.yGridChebyshev;
        end
    end
end