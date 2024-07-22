% © 2022 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
function phaseShiftVector = generate_phase_shift_vector(wavenumberObject, spatialShift)
    wavenumberVector = wavenumberObject.get_radial_frequency_vector();
    phaseShiftVector = exp(1.0j * wavenumberVector * spatialShift);
    nWavenumbers = numel(wavenumberVector);
    if mod(nWavenumbers, 2) == 0  % Nyquist frequency requires special treatment
        idxNyquistFreq = (nWavenumbers / 2) + 1;
        phaseShiftVector(idxNyquistFreq) = cos(wavenumberVector(idxNyquistFreq) * spatialShift);
    end
end