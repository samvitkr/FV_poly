% © 2023 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
function confTensorFourInterp = interpolate_confirmation_tensor_fourier_in_z(confTensorFour, wavenumberObjectZ, shiftZ)
    confTensorFourInterp.Cxx = interpolate_fourier_in_z(confTensorFour.Cxx, wavenumberObjectZ, shiftZ);
    confTensorFourInterp.Cyy = interpolate_fourier_in_z(confTensorFour.Cyy, wavenumberObjectZ, shiftZ);
    confTensorFourInterp.Czz = interpolate_fourier_in_z(confTensorFour.Czz, wavenumberObjectZ, shiftZ);
    confTensorFourInterp.Cxy = interpolate_fourier_in_z(confTensorFour.Cxy, wavenumberObjectZ, shiftZ);
    confTensorFourInterp.Cxz = interpolate_fourier_in_z(confTensorFour.Cxz, wavenumberObjectZ, shiftZ);
    confTensorFourInterp.Cyz = interpolate_fourier_in_z(confTensorFour.Cyz, wavenumberObjectZ, shiftZ);
end
