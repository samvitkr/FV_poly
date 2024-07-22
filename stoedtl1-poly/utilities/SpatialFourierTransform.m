% © 2022 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
classdef SpatialFourierTransform
    methods
        function [qFourier, wavenumberVectorX, wavenumberVectorZ] = transform_to_fourier(obj, qPhysical, dxSampling, dzSampling)
            arguments
                obj SpatialFourierTransform
                qPhysical double {mustBeReal}
                dxSampling double {mustBeReal} = -1.0  % default value ensures no wavenumber vector can be generated if argument is omitted
                dzSampling double {mustBeReal} = -1.0
            end
            % fourier transform
            if isvector(qPhysical)  % qPhysical is a vector of size 1xn or nx1: 1D transform
                nx = numel(qPhysical);
                normalizationFactor = 1.0 / nx;
                qFourier = normalizationFactor * fft(qPhysical);
            else
                nDims = ndims(qPhysical);
                if nDims == 2  % qPhysical is a matrix of size mxn (with m, n not equal to 1): 2D transform
                    [nz, nx] = size(qPhysical);
                    normalizationFactor = 1.0 / (nx * nz);
                    qFourier = normalizationFactor * fft2(qPhysical);
                elseif nDims == 3  % 3D data: repeated 2D transforms at each y location
                    [nz, nx, ny] = size(qPhysical);
                    normalizationFactor = 1.0 / (nx * nz);
                    qFourier = complex(zeros(nz, nx, ny));
                    for j = 1 : ny
                        qFourier(:, :, j) = normalizationFactor * fft2(qPhysical(:, :, j));
                    end
                else
                    error('fft not implemented for size of input array.');
                end
            end
            % generate wavenumber vectors (if requested)
            if nargout > 1  % in x
                if dxSampling >= 0.0
                    wavenumberVectorX = FrequencyVector(nx, dxSampling);
                else
                    error('invalid sampling rate in x.')
                end
            end
            if nargout > 2  % in z
                if dzSampling >= 0.0
                    wavenumberVectorZ = FrequencyVector(nz, dzSampling);
                else
                    error('invalid sampling rate in z.');
                end
            end
        end
        
        function [confirmationTensorFourier, wavenumberVectorX, wavenumberVectorZ] = transform_confirmation_tensor_to_fourier(obj, ...
                confirmationTensorPhysical, dxSampling, dzSampling)
            arguments
                obj SpatialFourierTransform
                confirmationTensorPhysical struct
                dxSampling double {mustBeReal} = -1.0  % default value ensures no wavenumber vector can be generated if argument is omitted
                dzSampling double {mustBeReal} = -1.0
            end
            % the wavenumber vectors are the same for all components, so compute them only once (if desired)
            if nargout == 1
                confirmationTensorFourier.Cxx = obj.transform_to_fourier(confirmationTensorPhysical.Cxx);
            elseif nargout == 2
                [confirmationTensorFourier.Cxx, wavenumberVectorX] = obj.transform_to_fourier(confirmationTensorPhysical.Cxx, ...
                    dxSampling);
            elseif nargout == 3
                [confirmationTensorFourier.Cxx, wavenumberVectorX, wavenumberVectorZ] = ...
                    obj.transform_to_fourier(confirmationTensorPhysical.Cxx, dxSampling, dzSampling);
            else
                error('Undefined number of output arguments');
            end
            % for subsequent components, we do not need to obtain the wavenumber vectors
            confirmationTensorFourier.Cyy = obj.transform_to_fourier(confirmationTensorPhysical.Cyy);
            confirmationTensorFourier.Czz = obj.transform_to_fourier(confirmationTensorPhysical.Czz);
            confirmationTensorFourier.Cxy = obj.transform_to_fourier(confirmationTensorPhysical.Cxy);
            confirmationTensorFourier.Cxz = obj.transform_to_fourier(confirmationTensorPhysical.Cxz);
            confirmationTensorFourier.Cyz = obj.transform_to_fourier(confirmationTensorPhysical.Cyz);
        end
        
        function qPhysical = transform_to_physical(obj, qFourier)
            arguments
                obj SpatialFourierTransform
                qFourier {mustBeNumeric}  % is there a way to enforce complex input?
            end
            if isvector(qFourier)  % qFourier is a vector of size 1xn or nx1: 1D transform
                nx = numel(qFourier);
                normalizationFactor = nx;
                qPhysical = normalizationFactor * ifft(qFourier, 'symmetric');
            else
                nDims = ndims(qFourier);
                if nDims == 2  % qFourier is a matrix of size mxn (with m, n not equal to 1): 2D transform
                    [nz, nx] = size(qFourier);
                    normalizationFactor = nx * nz;
                    qPhysical = normalizationFactor * ifft2(qFourier, 'symmetric');
                elseif nDims == 3  % 3D data: repeated 2D transforms at each y location
                    [nz, nx, ny] = size(qFourier);
                    normalizationFactor = nx * nz;
                    qPhysical = zeros(nz, nx, ny);
                    for j = 1 : ny
                        qPhysical(:, :, j) = normalizationFactor * ifft2(qFourier(:, :, j), 'symmetric');
                    end
                else
                    error('ifft not implemented for size of input array.');
                end     
            end
        end
        
        function confirmationTensorPhysical = transform_confirmation_tensor_to_physical(obj, confirmationTensorFourier)
            confirmationTensorPhysical.Cxx = obj.transform_to_physical(confirmationTensorFourier.Cxx);
            confirmationTensorPhysical.Cyy = obj.transform_to_physical(confirmationTensorFourier.Cyy);
            confirmationTensorPhysical.Czz = obj.transform_to_physical(confirmationTensorFourier.Czz);
            confirmationTensorPhysical.Cxy = obj.transform_to_physical(confirmationTensorFourier.Cxy);
            confirmationTensorPhysical.Cxz = obj.transform_to_physical(confirmationTensorFourier.Cxz);
            confirmationTensorPhysical.Cyz = obj.transform_to_physical(confirmationTensorFourier.Cyz);
        end
    end
end
