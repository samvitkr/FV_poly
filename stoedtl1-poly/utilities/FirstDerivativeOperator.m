% © 2022 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
classdef FirstDerivativeOperator
    properties (Access = protected)
        derivativeOperator {mustBeNumeric}
    end
    
    methods (Abstract)
        compute_derivative(obj, dataArray)
    end
    
    methods
        function derivativeOperator = get_derivative_operator(obj)
            derivativeOperator = obj.derivativeOperator;
        end
    end
end