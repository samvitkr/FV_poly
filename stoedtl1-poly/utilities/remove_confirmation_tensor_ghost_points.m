% © 2023 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
function confirmationTensorNoGhosts = remove_confirmation_tensor_ghost_points(confirmationTensorWithGhosts)
    confirmationTensorNoGhosts.Cxx = confirmationTensorWithGhosts.Cxx(2:end-1, 2:end-1, :);  % remove ghost points in x, z
    confirmationTensorNoGhosts.Cyy = confirmationTensorWithGhosts.Cyy(2:end-1, 2:end-1, :);
    confirmationTensorNoGhosts.Czz = confirmationTensorWithGhosts.Czz(2:end-1, 2:end-1, :);
    confirmationTensorNoGhosts.Cxy = confirmationTensorWithGhosts.Cxy(2:end-1, 2:end-1, :);
    confirmationTensorNoGhosts.Cxz = confirmationTensorWithGhosts.Cxz(2:end-1, 2:end-1, :);
    confirmationTensorNoGhosts.Cyz = confirmationTensorWithGhosts.Cyz(2:end-1, 2:end-1, :);
end
