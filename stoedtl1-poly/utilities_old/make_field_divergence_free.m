% © 2022 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
function [qXFourierDivFree, qYFourierDivFree, qZFourierDivFree] = make_field_divergence_free(qXFourier, qYFourier, ...
    qZFourier, wavenumberObjectX, wavenumberObjectZ, nCheb)
    % compute divergence residual
    wavenumberVectorX = wavenumberObjectX.get_radial_frequency_vector();
    wavenumberVectorZ = wavenumberObjectZ.get_radial_frequency_vector();
    dx = FirstDerivativeXFourier(wavenumberVectorX);
    dy = FirstDerivativeYChebyshev(nCheb);
    dz = FirstDerivativeZFourier(wavenumberVectorZ);
    residual = dx.compute_derivative(qXFourier) + dy.compute_derivative(qYFourier) + dz.compute_derivative(qZFourier);
    % allocate operators and arrays
    [~, DM] = chebdif(nCheb, 2);
    D1 = DM(:, :, 1);
    D2 = DM(:, :, 2);
    [nz, nx, ny] = size(qXFourier);
    scalarPotential = complex(zeros(nz, nx, ny));
    % transpose arrays for faster data access
    idxRearranged = [3 2 1];  % rearranges to y, x, z
    residual = permute(residual, idxRearranged);
    scalarPotential = permute(scalarPotential, idxRearranged);
    % solve Poisson equation and transpose result back
    for j = 1 : nz
        kz = wavenumberVectorZ(j);
        for i = 1 : nx
            if (j ~= 1 || i ~= 1)  % the Poisson problem for i = j = 1 is ill-posed. This mode is treated analytically below.
                kx = wavenumberVectorX(i);
                rhs = squeeze(residual(:, i, j));
                linearOperator = D2 - (kx^2 + kz^2) * eye(ny);
                % apply boundary conditions
                linearOperator(1, :) = D1(1, :);
                linearOperator(end, :) = D1(end, :);
                rhs(1) = complex(0.0);
                rhs(end) = complex(0.0);
                scalarPotential(:, i, j) = linearOperator \ rhs;
            end
        end
    end
    scalarPotential = ipermute(scalarPotential, idxRearranged);
    % compute irrotational part
    qXFourierIrrot = dx.compute_derivative(scalarPotential);
    qYFourierIrrot = dy.compute_derivative(scalarPotential);
    qZFourierIrrot = dz.compute_derivative(scalarPotential);
    % divergence-free part is the difference
    qXFourierDivFree = qXFourier - qXFourierIrrot;
    qYFourierDivFree = qYFourier - qYFourierIrrot;
    qZFourierDivFree = qZFourier - qZFourierIrrot;
    % treatment of the wall-parallel mean (i = j = 1, or, equivalently, kx = kz = 0), which was omitted above:
    % - x and z component: the mean x and z component of the irrotational part is zero by definition, so that the mean
    %                      of the solenoidal part coincides with the mean of the original data. This is correctly
    %                      enforced above, since scalarPotential(1, 1, :) == 0.
    % - y component: functional orthogonality (which ensures existence and uniqueness) requires that the y-component
    %                of the irroational part coincides with the y-component of the original data on the wall. The
    %                solenoidal part is therefore zero on the wall, which implies that its mean is zero on the wall also.
    %                The solenoidal constraint and periodic boundary conditions require that the mean of the solenoidal
    %                part is constant throughout the domain and coincides with the value at the wall, which is zero.
    %                We therefore set the y-component of the kx = kz = 0 mode explicity to zero below.
    qYFourierDivFree(1, 1, :) = 0.0;
end
