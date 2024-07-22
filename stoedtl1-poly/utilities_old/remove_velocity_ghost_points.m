% © 2022 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
function [uNoGhosts, vNoGhosts, wNoGhosts] = remove_velocity_ghost_points(uWithGhosts, vWithGhosts, wWithGhosts)
    uNoGhosts = uWithGhosts(2:end-1, 2:end-2, :);  % u on y cell centers. We need grid points 2:end-2 in x to make the data range consistent with v, w
    vNoGhosts = vWithGhosts(2:end-1, 2:end-1, 2:end-1);  % v on y cell vertices --> remove ghost points in y
    wNoGhosts = wWithGhosts(2:end-2, 2:end-1, :);  % w on y cell enters. The last gridpoint in z is a copy of the first, therefore 2:end-2
end