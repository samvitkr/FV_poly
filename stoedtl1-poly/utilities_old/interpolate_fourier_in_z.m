% © 2022 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
function qFourierInterpolatedZ = interpolate_fourier_in_z(qFourier, wavenumberObjectZ, shiftZ)
    phaseShift = generate_phase_shift_vector(wavenumberObjectZ, shiftZ);
    if isvector(qFourier)  % qFourier is a nx1 or 1xn vector
        % we need to make sure the first dimension of phase shift and qFourier are non-singleton. If the singleton
        % dimensions were mismatched, the multiplication below would result in an incorrect outer product.
        qFourierInterpolatedZ = shiftdim(phaseShift) .* shiftdim(qFourier);
    else
        nDims = ndims(qFourier);
        if nDims == 2  % qFourier is a mxn matrix (with m, n not equal to 1): 2d data
            [~, nx] = size(qFourier);
            phaseShift = repmat(reshape(phaseShift, [], 1), 1, nx);
            qFourierInterpolatedZ = phaseShift .* qFourier;
        elseif nDims == 3  % 3d data: apply phase shift at all y locations
            [~, nx, ny] = size(qFourier);
            phaseShift = repmat(reshape(phaseShift, [], 1, 1), 1, nx, ny);
            qFourierInterpolatedZ = phaseShift .* qFourier;
        else
            error('z interpolation not implemented for size of data array');
        end
    end
end