function dNdu = deriv_p(ipx, nnodel)
% Derivatives of shape functions with respect to local coordinates for one point
% $Id$

% local coordinates
xi  = ipx(1);
eta = ipx(2);

switch nnodel

    case 1

        % Derivatives of shape functions
        dNdu = [0.0; 0.0];

    case 3

        % Derivatives of shape functions
        dNdu = [
            % wrt xi
            0.0 1.0 0.0; ...
            % wrt eta
            0.0 0.0 1.0];

    case 4

        % 1D Polynomials
        p10 = 0.5 * (1.0 - xi);
        p11 = 0.5 * (1.0 + xi);
        p20 = 0.5 * (1.0 - eta);
        p21 = 0.5 * (1.0 + eta);

        % Derivatives of 1D Polynomials
        p10d = -0.5;
        p11d =  0.5;
        p20d = -0.5;
        p21d =  0.5;

        % Derivatives of shape functions
        dNdu = [
            %wrt xi
            p10d * p20 ...
            p11d * p20 ...
            p11d * p21 ...
            p10d * p21; ...
            %wrt eta
            p10 * p20d ...
            p11 * p20d ...
            p11 * p21d ...
            p10 * p21d];

    case 9

        % 1D Polynomials
        p10 = 0.5 * xi * (xi - 1.0);
        p11 = -(xi + 1.0) * (xi - 1.0);
        p12 = 0.5 * xi * (xi + 1.0);
        p20 = 0.5 * eta * (eta - 1.0);
        p21 = -(eta + 1.0) * (eta - 1.0);
        p22 = 0.5 * eta * (eta + 1.0);

        % Derivatives of 1D Polynomials
        p10d = 0.5 * (2 * xi - 1.0);
        p11d = -2.0 * xi;
        p12d = 0.5 * (2 * xi + 1.0);
        p20d = 0.5 * (2 * eta - 1.0);
        p21d = -2.0 * eta;
        p22d = 0.5 * (2 * eta + 1.0);

        % Derivatives of shape functions
        dNdu = [
            %wrt xi
            p10d * p20 ...
            p12d * p20 ...
            p12d * p22 ...
            p10d * p22 ...
            p11d * p20 ...
            p12d * p21 ...
            p11d * p22 ...
            p10d * p21 ...
            p11d * p21; ...
            %wrt eta
            p10 * p20d ...
            p12 * p20d ...
            p12 * p22d ...
            p10 * p22d ...
            p11 * p20d ...
            p12 * p21d ...
            p11 * p22d ...
            p10 * p21d ...
            p11 * p21d];

    otherwise
        error('Unknown element')

end

end