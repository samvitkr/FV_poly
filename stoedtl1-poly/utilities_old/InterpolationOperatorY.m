% © 2022 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
classdef InterpolationOperatorY
    properties (Access = protected)
        % note: the y-coordinate range is taken as [-1, 1], as required by the Chebyshev grid.
        % The constructor enforces this y range explicitly.
        yGridPointsCheb;
        yGridPointsDns;
        yCellCentersDns;
    end
    
    methods
        function obj = InterpolationOperatorY(yGridPointsDns, yGridPointsCheb)
            yGridPointsDns = shiftdim(yGridPointsDns);  % make sure first dimension is non-singleton
            yGridPointsCheb = shiftdim(yGridPointsCheb);
            absTol = 1e-8;
            % make sure the Chebyshev grid is [-1, 1], in descending order
            differenceTop = yGridPointsCheb(1) - 1.0;
            differenceBot = yGridPointsCheb(end) + 1.0;
            if (abs(differenceBot) > absTol) || (abs(differenceTop) > absTol)
                error('Y interpolation not implemented for given y Chebyshev grid')
            end        
            obj.yGridPointsCheb = yGridPointsCheb;
            % shift the DNS grid if necessary to make its range [-1, 1]
            differenceToChebBot = yGridPointsDns(1) + 1.0;
            differenceToChebTop = yGridPointsDns(end) - 1.0;
            differenceToFiniteVolBot = yGridPointsDns(1);
            differenceToFiniteVolTop = yGridPointsDns(end) - 2.0;
            if (abs(differenceToChebBot) > absTol) || (abs(differenceToChebTop) > absTol)  % DNS grid is not [-1, 1]
                if (abs(differenceToFiniteVolBot) < absTol) && (abs(differenceToFiniteVolTop) < absTol)  % DNS grid is [0, 2]
                    yGridPointsDns = yGridPointsDns - 1.0;  % shift to [-1, 1]
                else
                    error('Y interpolation not implemented for range of y grid.');
                end
            end
            obj.yGridPointsDns = yGridPointsDns;
            yCellCentersDns = 0.5 * (yGridPointsDns(1:end-1) + yGridPointsDns(2:end));  % interior points
            obj.yCellCentersDns = [yGridPointsDns(1); yCellCentersDns; yGridPointsDns(end)];  % add wall location
        end
        
        function uCheb = interpolate_u_to_chebyshev_grid(obj, uCenterWithGhostPoints)
            uCenterNoSlip = obj.replace_ghost_points_with_no_slip_bc(uCenterWithGhostPoints);
            uCheb = obj.interpolate_from_cell_centers_to_chebyshev_points(uCenterNoSlip);
        end
        
        function vCheb = interpolate_v_to_chebyshev_grid(obj, vGridPoints)
            vWithSlopeConstraint = obj.add_slope_constraint(vGridPoints);  % dv/dy = 0 at the wall
            vCheb = obj.interpolate_from_grid_points_to_chebyshev_points(vWithSlopeConstraint);
        end
        
        function wCheb = interpolate_w_to_chebyshev_grid(obj, wCenterWithGhostPoints)
            wCenterNoSlip = obj.replace_ghost_points_with_no_slip_bc(wCenterWithGhostPoints);
            wCheb = obj.interpolate_from_cell_centers_to_chebyshev_points(wCenterNoSlip);
        end
        
        function TCheb = interpolate_temperature_to_chebyshev_grid(obj, TCenterWithGhostPoints, TBottom, TTop)
            arguments
                obj InterpolationOperatorY
                TCenterWithGhostPoints {mustBeNumeric}  % input can be real or complex
                TBottom {mustBeNumeric} = 0.0
                TTop {mustBeNumeric} = 0.0
            end
            TCenterWithBc = obj.replace_ghost_points_with_dirichlet_bc(TCenterWithGhostPoints, TBottom, TTop);
            TCheb = obj.interpolate_from_cell_centers_to_chebyshev_points(TCenterWithBc);
        end
        
        function confTensorCheb = interpolate_confirmation_tensor_to_chebyshev_grid(obj, confTensorWithGhostPoints)
            confTensorCheb.Cxx = obj.interpolate_tensor_component_to_chebyshev_points(confTensorWithGhostPoints.Cxx);
            confTensorCheb.Cyy = obj.interpolate_tensor_component_to_chebyshev_points(confTensorWithGhostPoints.Cyy);
            confTensorCheb.Czz = obj.interpolate_tensor_component_to_chebyshev_points(confTensorWithGhostPoints.Czz);
            confTensorCheb.Cxy = obj.interpolate_tensor_component_to_chebyshev_points(confTensorWithGhostPoints.Cxy);
            confTensorCheb.Cxz = obj.interpolate_tensor_component_to_chebyshev_points(confTensorWithGhostPoints.Cxz);
            confTensorCheb.Cyz = obj.interpolate_tensor_component_to_chebyshev_points(confTensorWithGhostPoints.Cyz);
        end
        
        function yGridPointsCheb = get_y_grid_points_chebyshev(obj)
            yGridPointsCheb = obj.yGridPointsCheb;
        end
        
        function yGridPointsDns = get_y_grid_points_dns(obj)
            yGridPointsDns = obj.yGridPointsDns;
        end
        
        function yCellCentersDns = get_y_cell_centers_dns(obj)
            yCellCentersDns = obj.yCellCentersDns;
        end
    end
    
    methods (Access = protected)
        function qWithBc = replace_ghost_points_with_no_slip_bc(obj, qWithGhostPoints)
            qWithBc = obj.replace_ghost_points_with_dirichlet_bc(qWithGhostPoints, 0.0, 0.0);
        end
        
        function qWithBc = replace_ghost_points_with_dirichlet_bc(obj, qWithGhostPoints, dirichletDataBot, dirichletDataTop)
            % note: the y ordering of the data matters when setting general Dirichlet boundary conditions. qWithBc is
            % assumed to live on the DNS grid, so that y data order is [-1 -> +1], i.e. in increasing order
            qWithBc = qWithGhostPoints;  % copy data over since interior data will be unchanged
            if isvector(qWithBc)  % qWithBc is 1xn or nx1 vector
                qWithBc(1) = dirichletDataBot;  % bottom wall
                qWithBc(end) = dirichletDataTop;  % top wall
            else
                nDims = ndims(qWithBc);
                if nDims == 3
                    qWithBc(:, :, 1) = dirichletDataBot;  % bottom wall
                    qWithBc(:, :, end) = dirichletDataTop;  % top wall
                else
                    error('assignment of Dirichlet boundary condition not implemented for size of data array');
                end
            end
        end
        
        function qWithSlopeConstraint = add_slope_constraint(obj, qWithBc)
            if isvector(qWithBc)
                n = numel(qWithBc);
                qWithSlopeConstraint = zeros(n+2, 1);
                qWithSlopeConstraint(2:end-1) = qWithBc;
            else
                nDims = ndims(qWithBc);
                if nDims == 3
                    [nz, nx, ny] = size(qWithBc);
                    qWithSlopeConstraint = zeros(nz, nx, ny+2);
                    qWithSlopeConstraint(:, :, 2:end-1) = qWithBc;
                else
                    error('assignment of slope constraint not implemented for size of data array');
                end
            end
        end
        
        function qChebyshev = interpolate_from_cell_centers_to_chebyshev_points(obj, qCenter)
            % make sure qCenter is a vector or 3D array and retreive y-dimension
            if isvector(qCenter)
                ny = numel(qCenter);
            else
                nDims = ndims(qCenter);
                if nDims == 3
                    [~, ~, ny] = size(qCenter);
                else  % note: unlike for x and z, the notion of 2D data is not meaningful for the y-coordinate
                    error('Y interpolation not implemented for size of input array');
                end
            end
            % only do interpolation if qCenter has the correct y-dimensions (possibly including a slope constraint)
            if ny == numel(obj.yCellCentersDns) || ny == (numel(obj.yCellCentersDns) + 2)
                % interpolate in y. spline works on the last dimension of qCenter, which always corresponds to the
                % y-coordinate. We can therefore treat the 1D and 3D case together.
                % Note: qChebyshev inherits the data ordering of yGridPointsCheb, i.e. descending in y
                qChebyshev = spline(obj.yCellCentersDns, qCenter, obj.yGridPointsCheb);
            else
                error('Y interpolation not implemented for size of input array');
            end
        end
        
        function qChebyshev = interpolate_from_grid_points_to_chebyshev_points(obj, qGridPoints)
            % make sure qGridPoint is a vector or 3D array and retreive y-dimension
            if isvector(qGridPoints)
                ny = numel(qGridPoints);
            else
                nDims = ndims(qGridPoints);
                if nDims == 3
                    [~, ~, ny] = size(qGridPoints);
                else
                    error('Y interpolation not implemented for size of input array');
                end
            end
            % only do interpolation if qGridPoints has the correct y-dimensions (possibly including slope constraint)
            if ny == (numel(obj.yGridPointsDns) + 2) || ny == numel(obj.yGridPointsDns)
                % interpolate in y. spline works on the last dimension of qCenter, which always corresponds to the
                % y-coordinate. We can therefore treat the 1D and 3D case together.
                % Note: qChebyshev inherits the data ordering of yGridPointsCheb, i.e. descending in y
                qChebyshev = spline(obj.yGridPointsDns, qGridPoints, obj.yGridPointsCheb);
            else
                error('Y interpolation not implemented for size of input array');
            end
        end
        
        function tensorComponentCheb = interpolate_tensor_component_to_chebyshev_points(obj, tensorComponentWithGhostPoints)
            % interpolate wall value of confirmation tensor component and set boundary conditions
            if isvector(tensorComponentWithGhostPoints)
                dirichletDataBot = 0.5 * (tensorComponentWithGhostPoints(1) + tensorComponentWithGhostPoints(2));
                dirichletDataTop = 0.5 * (tensorComponentWithGhostPoints(end-1) + tensorComponentWithGhostPoints(end));
            else
                nDims = ndims(tensorComponentWithGhostPoints);
                if nDims == 3
                    dirichletDataBot = 0.5 * (tensorComponentWithGhostPoints(:, :, 1) ...
                        + tensorComponentWithGhostPoints(:, :, 2));
                    dirichletDataTop = 0.5 * (tensorComponentWithGhostPoints(:, :, end-1) ...
                        + tensorComponentWithGhostPoints(:, :, end));
                else
                    error('Y interpolation not implemented for size of input array');
                end
            end
            tensorComponentWithBc = obj.replace_ghost_points_with_dirichlet_bc(tensorComponentWithGhostPoints, ...
                dirichletDataBot, dirichletDataTop);
            tensorComponentCheb = obj.interpolate_from_cell_centers_to_chebyshev_points(tensorComponentWithBc);
        end
    end
end