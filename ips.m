function [ipx, ipw] = ips(nip)
% Integration rules (points & weights)
% $Id$

switch nip
    case 4
        % Local coordinates of integration points (xi,eta)
        g1 = -sqrt(1/3);
        g2 =  sqrt(1/3);
        ipx(1,1:2) = [g1, g1];
        ipx(2,1:2) = [g2, g1];
        ipx(3,1:2) = [g2, g2];
        ipx(4,1:2) = [g1, g2];

        % Corresponding weights
        w = 1.0;
        ipw(1) = w;
        ipw(2) = w;
        ipw(3) = w;
        ipw(4) = w;

    case 9
        % Local coordinates of integration points (xi,eta)
        g1 = -sqrt(3/5);
        g2 =  sqrt(3/5);
        ipx(1,1:2) = [g1, g1];
        ipx(2,1:2) = [g2, g1];
        ipx(3,1:2) = [g2, g2];
        ipx(4,1:2) = [g1, g2];
        ipx(5,1:2) = [0,  g1];
        ipx(6,1:2) = [g2,  0];
        ipx(7,1:2) = [0,  g2];
        ipx(8,1:2) = [g1,  0];
        ipx(9,1:2) = [0,   0];

        % Corresponding weights
        w1 = 5/9;
        w2 = 8/9;
        w11 = w1 * w1;
        w12 = w1 * w2;
        w22 = w2 * w2;
        ipw(1) = w11;
        ipw(2) = w11;
        ipw(3) = w11;
        ipw(4) = w11;
        ipw(5) = w12;
        ipw(6) = w12;
        ipw(7) = w12;
        ipw(8) = w12;
        ipw(9) = w22;

    otherwise
        error('Unknown integration rule')

end

end