% © 2022 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
classdef FirstDerivativeOperatorFourier < FirstDerivativeOperator
    methods
        function obj = FirstDerivativeOperatorFourier(wavenumberVector)
            nWavenumbers = numel(wavenumberVector);
            obj.derivativeOperator = 1.0j * wavenumberVector;
            if mod(nWavenumbers, 2) == 0  % Nyquist frequency is present only if wavenumber vector has an even length
                idxNyquistFreq = (nWavenumbers / 2) + 1;
                % the first deriviative of the Nyuist frequency is zero on the discretized grid (important note: this is
                % not necessarily true for higher-order derivatives)
                obj.derivativeOperator(idxNyquistFreq) = complex(0);
            end
        end
    end
end