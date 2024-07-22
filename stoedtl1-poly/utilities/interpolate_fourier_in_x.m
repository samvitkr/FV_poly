% © 2022 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
function qFourierInterpolatedX = interpolate_fourier_in_x(qFourier, wavenumberObjectX, shiftX)
    phaseShift = generate_phase_shift_vector(wavenumberObjectX, shiftX);
    % apply phase shift
    if isvector(qFourier)  % qFourier is a nx1 or 1xn vector
        % we need to make sure the first dimension of phase shift and qFourier are non-singleton. If the singleton
        % dimensions were mismatched, the multiplication below would result in an incorrect outer product. 
        qFourierInterpolatedX = shiftdim(phaseShift) .* shiftdim(qFourier);
    else
        nDims = ndims(qFourier);
        if nDims == 2
            [nz, ~] = size(qFourier);
            phaseShift = repmat(reshape(phaseShift, 1, []), nz, 1);
            qFourierInterpolatedX = phaseShift .* qFourier;
        elseif nDims == 3
            [nz, ~, ny] = size(qFourier);
            phaseShift = repmat(reshape(phaseShift, 1, [], 1), nz, 1, ny);
            qFourierInterpolatedX = phaseShift .* qFourier;
        else
            error('x interpolation not implemented for size of data array.')
        end
    end
end