function N = shp_p(ipx, nnodel)
% Shape functions with respect to local coordinates for one point
% $Id$

% local coordinates
xi  = ipx(1);
eta = ipx(2);

switch nnodel

    case 1

        % Shape functions
        N = [1.0];

    case 3

        % Shape functions
        N = [1.0 xi eta];

    case 4

        % 1D Polynomials
        p10 = 0.5 * (1.0 - xi);
        p11 = 0.5 * (1.0 + xi);
        p20 = 0.5 * (1.0 - eta);
        p21 = 0.5 * (1.0 + eta);

        % Shape functions
        N = [
            p10 * p20 ...
            p11 * p20 ...
            p11 * p21 ...
            p10 * p21];

    case 9

        % 1D Polynomials
        p10 = 0.5 * xi * (xi - 1.0);
        p11 = -(xi + 1.0) * (xi - 1.0);
        p12 = 0.5 * xi * (xi + 1.0);
        p20 = 0.5 * eta * (eta - 1.0);
        p21 = -(eta + 1.0) * (eta - 1.0);
        p22 = 0.5 * eta * (eta + 1.0);

        % Shape functions
        N = [
            p10 * p20 ...
            p12 * p20 ...
            p12 * p22 ...
            p10 * p22 ...
            p11 * p20 ...
            p12 * p21 ...
            p11 * p22 ...
            p10 * p21 ...
            p11 * p21];

    otherwise
        error('Unknown element')

end

end