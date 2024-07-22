% © 2023 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
function confTensorFourInterp = interpolate_confirmation_tensor_fourier_in_x(confTensorFour, wavenumberObjectX, shiftX)
    confTensorFourInterp.Cxx = interpolate_fourier_in_x(confTensorFour.Cxx, wavenumberObjectX, shiftX);
    confTensorFourInterp.Cyy = interpolate_fourier_in_x(confTensorFour.Cyy, wavenumberObjectX, shiftX);
    confTensorFourInterp.Czz = interpolate_fourier_in_x(confTensorFour.Czz, wavenumberObjectX, shiftX);
    confTensorFourInterp.Cxy = interpolate_fourier_in_x(confTensorFour.Cxy, wavenumberObjectX, shiftX);
    confTensorFourInterp.Cxz = interpolate_fourier_in_x(confTensorFour.Cxz, wavenumberObjectX, shiftX);
    confTensorFourInterp.Cyz = interpolate_fourier_in_x(confTensorFour.Cyz, wavenumberObjectX, shiftX);
end
