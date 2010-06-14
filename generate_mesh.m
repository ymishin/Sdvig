function [X, g_num, bpoints] = generate_mesh(xyvec, el, nodes, no_nodes_el)
% $Id$

% number of elements
no_el_x = el(1);            % along x direction
no_el_y = el(2);            % along y direction
no_el = no_el_x * no_el_y;  % total

% number of nodes
no_nodes_x = nodes(1);              % along x direction
no_nodes_y = nodes(2);              % along y direction
no_nodes = no_nodes_x * no_nodes_y; % total

% global coordinates of nodes (up-right numbering)
xvec = xyvec{1};
yvec = xyvec{2};
xmat = reshape(repmat(xvec,no_nodes_y,1),1,no_nodes);
ymat = repmat(yvec,1,no_nodes_x);
X = [xmat; ymat];

% boundary points
left    = (1:no_nodes_y);
right   = (no_nodes-no_nodes_y+1:no_nodes);
top     = (no_nodes_y:no_nodes_y:no_nodes);
bottom  = (1:no_nodes_y:no_nodes);
bpoints = {left right top bottom};

% build mesh (assign nodes to elements)
g_num = zeros(no_nodes_el,no_el);
switch no_nodes_el    
    case 1
        shft = repmat([0:no_nodes_y:no_nodes-no_nodes_y],no_el_y,1);
        % node #1
        n = repmat([1:1:no_nodes_y]',1,no_el_x) + shft;
        g_num(1,:) = reshape(n,1,no_el);
        
    case 4
        shft1 = repmat([0:no_nodes_y:no_nodes-2*no_nodes_y],no_el_y,1);
        shft2 = repmat([no_nodes_y:no_nodes_y:no_nodes-no_nodes_y],no_el_y,1);
        % node #1
        n = repmat([1:1:no_nodes_y-1]',1,no_el_x) + shft1;
        g_num(1,:) = reshape(n,1,no_el);
        % node #2
        n = repmat([1:1:no_nodes_y-1]',1,no_el_x) + shft2;
        g_num(2,:) = reshape(n,1,no_el);
        % node #3
        n = repmat([2:1:no_nodes_y]',1,no_el_x) + shft2;
        g_num(3,:) = reshape(n,1,no_el);
        % node #4
        n = repmat([2:1:no_nodes_y]',1,no_el_x) + shft1;
        g_num(4,:) = reshape(n,1,no_el);

    case 9
        shft1 = repmat([0:2*no_nodes_y:no_nodes-3*no_nodes_y],no_el_y,1);
        shft2 = repmat([2*no_nodes_y:2*no_nodes_y:no_nodes-no_nodes_y],no_el_y,1);
        shft3 = repmat([no_nodes_y:2*no_nodes_y:no_nodes-2*no_nodes_y],no_el_y,1);
        % node #1
        n = repmat([1:2:no_nodes_y-2]',1,no_el_x) + shft1;
        g_num(1,:) = reshape(n,1,no_el);
        % node #2
        n = repmat([1:2:no_nodes_y-2]',1,no_el_x) + shft2;
        g_num(2,:) = reshape(n,1,no_el);
        % node #3
        n = repmat([3:2:no_nodes_y]',1,no_el_x) + shft2;
        g_num(3,:) = reshape(n,1,no_el);
        % node #4
        n = repmat([3:2:no_nodes_y]',1,no_el_x) + shft1;
        g_num(4,:) = reshape(n,1,no_el);
        % node #5
        n = repmat([1:2:no_nodes_y-2]',1,no_el_x) + shft3;
        g_num(5,:) = reshape(n,1,no_el);
        % node #6
        n = repmat([2:2:no_nodes_y-1]',1,no_el_x) + shft2;
        g_num(6,:) = reshape(n,1,no_el);
        % node #7
        n = repmat([3:2:no_nodes_y]',1,no_el_x) + shft3;
        g_num(7,:) = reshape(n,1,no_el);
        % node #8
        n = repmat([2:2:no_nodes_y-1]',1,no_el_x) + shft1;
        g_num(8,:) = reshape(n,1,no_el);
        % node #9
        n = repmat([2:2:no_nodes_y-1]',1,no_el_x) + shft3;
        g_num(9,:) = reshape(n,1,no_el);

    otherwise
        error('Unknown mesh')

end

g_num = uint32(g_num);

end
