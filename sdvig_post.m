function sdvig_post
% $Id$

close all;
clear all;
clc;

% case name and counters
case_name = 'extension';
ifirst = 1;
ilast  = 1;
step = 1;

% material field (from particles)
plot_materialf = true;
save_materialf = false;
% strain2 field (from particles)
plot_strain2f  = true;
save_strain2f  = false;
% pressure field (from grid)
plot_pressuref = true;
save_pressuref = false;
% x- and y-velocity fields (from grid)
plot_velocityf = true;
save_velocityf = false;
plot_all_vnodes = false; % all velocity nodes OR only elements' corners
plot_vgrid = false;
% just velocity grid
plot_vgrid_only = false;

for i = ifirst:step:ilast
    
    % file name
    fname = [case_name,'_',num2str(i,'%04d'),'.h5'];
    disp(fname);
    
    % read data from the file
    time = hdf5read(fname,'/Info/CurrentTime');
    domain = hdf5read(fname,'/Info/InitialDomain');
    xmin0 = domain(1);
    xmax0 = domain(2);
    ymin0 = domain(3);
    ymax0 = domain(4);
    domain = hdf5read(fname,'/Info/Domain');
    xmin = domain(1);
    xmax = domain(2);
    ymin = domain(3);
    ymax = domain(4);
    res = hdf5read(fname,'/Info/InitialResolution');
    no_el_x = res(1);
    no_el_y = res(2);
    no_el = no_el_x * no_el_y;
    el_type = hdf5read(fname,'/Info/ElementType');
    part_per_side = hdf5read(fname,'/Info/ParticlesPerSide');
    part_xcoord = hdf5read(fname,'/Particles/XCoordinates');
    part_ycoord = hdf5read(fname,'/Particles/YCoordinates');
    part_material = hdf5read(fname,'/Particles/Material');
    part_strain2 = hdf5read(fname,'/Particles/StrainInvariant2');
    grid_vnodes = hdf5read(fname,'/Grid/VelocityNodes');
    grid_vnodes2elements = hdf5read(fname,'/Grid/VelocityNodes2Elements');
    grid_xvelocity = hdf5read(fname,'/Grid/XVelocity');
    grid_yvelocity = hdf5read(fname,'/Grid/YVelocity');
    grid_pressure = hdf5read(fname,'/Grid/Pressure');
    
    switch el_type
        case 1
            % Q1P0
            pressure = grid_pressure;
        case 2
            % Q1Q1
            pressure = reshape(grid_pressure(grid_vnodes2elements(1:4,:)), ...
                4,no_el);
            pressure = pressure';
        case 3
            % Q2P-1
            pressure = grid_pressure(1:3:end);
            pressure_gx = grid_pressure(2:3:end); % X gradient
            pressure_gy = grid_pressure(3:3:end); % Y gradient
    end
    
    % number of cells per element for interpolation
    cells_element = part_per_side + 1;
    
    % assign particles to grid cells (cells are numbered up-right)
    no_xcells = no_el_x * cells_element;
    no_ycells = no_el_y * cells_element;
    no_cells = no_xcells * no_ycells;
    cell_xsize = (xmax - xmin) / no_xcells;
    cell_xsize1 = 1.0 / cell_xsize;
    cell_ysize = (ymax - ymin) / no_ycells;
    cell_ysize1 = 1.0 / cell_ysize;
    part2cells = cell(no_ycells,no_xcells);
    no_part = size(part_xcoord,1);
    for ip = 1:no_part
        % check if particle is inside the domain
        if (part_xcoord(ip) < xmin || part_xcoord(ip) > xmax || ...
                part_ycoord(ip) < ymin || part_ycoord(ip) > ymax)
            continue;
        end
        % find the cell particle resides in
        k = ceil((part_xcoord(ip) - xmin) * cell_xsize1);
        l = ceil((part_ycoord(ip) - ymin) * cell_ysize1);
        part2cells{l,k} = [part2cells{l,k} ip];
    end
    
    % coordinates of the cells (TODO: remove - call external mesh generator)
    cgrid_xnodes = zeros(4,no_cells);
    cgrid_ynodes = zeros(4,no_cells);
    vxvec = zeros(no_xcells + 1,1);
    vyvec = zeros(no_ycells + 1,1);
    centers_xcells = zeros(no_xcells,1);
    centers_ycells = zeros(no_ycells,1);
    for k = 1:no_xcells
        centers_xcells(k) = xmin + (k - 0.5) * cell_xsize;
        vxvec(k) = xmin + (k - 1) * cell_xsize;
    end
    vxvec(k+1) = xmax;
    for k = 1:no_ycells
        centers_ycells(k) = ymin + (k - 0.5) * cell_ysize;
        vyvec(k) = ymin + (k - 1) * cell_ysize;
    end
    vyvec(k+1) = ymax;
    
    % fields to be interpolated from particles
    emptycells_field = zeros(no_ycells,no_xcells);
    material_field = zeros(no_ycells,no_xcells);
    strain2_field = zeros(no_ycells,no_xcells);
    
    % exclude empty area (from the top)
    for k = 1:no_xcells
        for l = no_ycells:-1:1
            % particles in the current cell
            parts = part2cells{l,k};
            if (~isempty(parts))
                break;
            end
            % mark "empty" cells
            emptycells_field(l,k) = 1;
        end
    end
    for k = 1:no_xcells
        for l = 1:no_ycells
            % ignore narrow lines of empty cells
            if (emptycells_field(l,k) && ...
                    ((k == 1) || ~emptycells_field(l,k-1)) && ...
                    ((k == no_xcells) || ~emptycells_field(l,k+1)))
                emptycells_field(l,k) = 0;
            end
        end
    end
    
    % interpolation from particles
    for k = 1:no_xcells
        for l = 1:no_ycells
            % coordinates of the cells (TODO: remove - call external mesh generator)
            cell_num = (k - 1) * no_ycells + l;
            cgrid_xnodes(1,cell_num) = vxvec(k);
            cgrid_xnodes(2,cell_num) = vxvec(k+1);
            cgrid_xnodes(3,cell_num) = vxvec(k+1);
            cgrid_xnodes(4,cell_num) = vxvec(k);
            cgrid_ynodes(1,cell_num) = vyvec(l);
            cgrid_ynodes(2,cell_num) = vyvec(l);
            cgrid_ynodes(3,cell_num) = vyvec(l+1);
            cgrid_ynodes(4,cell_num) = vyvec(l+1);
            % check if the current cell is "empty"
            if (emptycells_field(l,k))
                continue;
            end
            % particles in the current cell
            parts = part2cells{l,k};
            % is the current cell a hole?
            if (isempty(parts))
                m = 1;
                % look at the neighbour cells
                while (true)
                    for k1 = -m:1:m
                        % check domain boundaries
                        if ((k + k1) < 1 || (k + k1) > no_xcells)
                            continue;
                        end
                        if (k1 == -m || k1 == m)
                            ld = [-m:1:m];
                        else
                            ld = [-m m];
                        end
                        for l1 = ld
                            % check domain boundaries
                            if ((l + l1) < 1 || (l + l1) > no_ycells)
                                continue;
                            end
                            parts = [parts part2cells{l+l1,k+k1}];
                        end
                    end
                    if (~isempty(parts))
                        break;
                    end
                    m = m + 1;
                end
            end
            if (length(parts) > 1)
                % find nearest particle to the cell center
                d = realmax;
                for ip = parts
                    dn = sqrt((centers_xcells(k) - part_xcoord(ip))^2 + ...
                        (centers_ycells(l) - part_ycoord(ip))^2);
                    if (dn < d)
                        interp_p = ip;
                        d = dn;
                    end
                end
            else
                % only one particle
                interp_p  = parts(1);
            end
            % interpolate
            material_field(l,k) = part_material(interp_p);
            strain2_field(l,k) = part_strain2(interp_p);
        end
    end
    
    % plot and save
    
    % domain to plot in
    axis_domain = [xmin xmax ymin ymax];
    
    % non-empty cells
    emptycells_field = reshape(emptycells_field,no_cells,1);
    nonempty_cells = find(emptycells_field ~= 1);
    
    % material field (from particles)
    if (plot_materialf)
        
        figure;
        c_map = [  1.0  1.0  1.0; ...
                  0.35 0.35 0.35; ...
                  0.50 0.50 0.50];
        colormap(c_map);
        caxis([1 3]);
        material_field = reshape(material_field,no_cells,1);
        patch(cgrid_xnodes(:,nonempty_cells),cgrid_ynodes(:,nonempty_cells), ...
            material_field(nonempty_cells)','EdgeColor','none');
        title('Material');
        box on;
        axis xy equal;
        axis(axis_domain);
        set(gca,'Layer','top');
        
        % save figure
        if (save_materialf)
            fname = [case_name,'_material_',num2str(i,'%04d')];
            print('-djpeg','-r200',[fname,'.jpg']);
            close;
        end
        
    end
    
    % strain2 field (from particles)
    if (plot_strain2f)
        
        figure;
        strain2_field = reshape(strain2_field,no_cells,1);
        patch(cgrid_xnodes(:,nonempty_cells),cgrid_ynodes(:,nonempty_cells), ...
            log10(strain2_field(nonempty_cells))','EdgeColor','none');
        title('Strain2');
        box on;
        axis xy equal;
        axis(axis_domain);
        set(gca,'Layer','top');
        
        % save figure
        if (save_strain2f)
            fname = [case_name,'_strain2_',num2str(i,'%04d')];
            print('-djpeg','-r200',[fname,'.jpg']);
            close;
        end
    end
    
    % pressure field (from grid)
    if (plot_pressuref)
        
        figure;
        patch(reshape(grid_vnodes(1,grid_vnodes2elements(1:4,:)),4,no_el), ...
            reshape(grid_vnodes(2,grid_vnodes2elements(1:4,:)),4,no_el), ...
            pressure','EdgeColor','none');
        title('Pressure');
        box on;
        axis xy equal;
        axis(axis_domain);
        set(gca,'Layer','top');
        
        % save figure
        if (save_pressuref)
            fname = [case_name,'_pressure_',num2str(i,'%04d')];
            print('-djpeg','-r200',[fname,'.jpg']);
            close;
        end
    end
    
    % velocity fields (from grid)
    if (plot_velocityf)
        
        % plot velocity grid ?
        if (plot_vgrid)
            edge_color = 'k';
        else
            edge_color = 'none';
        end
        
        if (plot_all_vnodes && (el_type == 3))
            % consider all nodes; only for Q2P-1 element (9 velocity nodes)
            coord_x = [reshape(grid_vnodes(1,grid_vnodes2elements([1 5 9 8],:)),4,no_el) ...
                reshape(grid_vnodes(1,grid_vnodes2elements([5 2 6 9],:)),4,no_el) ...
                reshape(grid_vnodes(1,grid_vnodes2elements([9 6 3 7],:)),4,no_el) ...
                reshape(grid_vnodes(1,grid_vnodes2elements([8 9 7 4],:)),4,no_el)];
            coord_y = [reshape(grid_vnodes(2,grid_vnodes2elements([1 5 9 8],:)),4,no_el) ...
                reshape(grid_vnodes(2,grid_vnodes2elements([5 2 6 9],:)),4,no_el) ...
                reshape(grid_vnodes(2,grid_vnodes2elements([9 6 3 7],:)),4,no_el) ...
                reshape(grid_vnodes(2,grid_vnodes2elements([8 9 7 4],:)),4,no_el)];
            vel_x = [reshape(grid_xvelocity(grid_vnodes2elements([1 5 9 8],:)),4,no_el) ...
                reshape(grid_xvelocity(grid_vnodes2elements([5 2 6 9],:)),4,no_el) ...
                reshape(grid_xvelocity(grid_vnodes2elements([9 6 3 7],:)),4,no_el) ...
                reshape(grid_xvelocity(grid_vnodes2elements([8 9 7 4],:)),4,no_el)];
            vel_y = [reshape(grid_yvelocity(grid_vnodes2elements([1 5 9 8],:)),4,no_el) ...
                reshape(grid_yvelocity(grid_vnodes2elements([5 2 6 9],:)),4,no_el) ...
                reshape(grid_yvelocity(grid_vnodes2elements([9 6 3 7],:)),4,no_el) ...
                reshape(grid_yvelocity(grid_vnodes2elements([8 9 7 4],:)),4,no_el)];
        else
            % consider nodes only at elements' corners
            coord_x = reshape(grid_vnodes(1,grid_vnodes2elements(1:4,:)),4,no_el);
            coord_y = reshape(grid_vnodes(2,grid_vnodes2elements(1:4,:)),4,no_el);
            vel_x = reshape(grid_xvelocity(grid_vnodes2elements(1:4,:)),4,no_el);
            vel_y = reshape(grid_yvelocity(grid_vnodes2elements(1:4,:)),4,no_el);
        end
        
        % x-velocity
        figure;
        patch(coord_x,coord_y,vel_x,'EdgeColor',edge_color);
        title('X-Velocity');
        box on;
        axis xy equal;
        axis(axis_domain);
        set(gca,'Layer','top');
        
        % save figure
        if (save_velocityf)
            fname = [case_name,'_x-velocity_',num2str(i,'%04d')];
            print('-djpeg','-r200',[fname,'.jpg']);
            close;
        end
        
        % y-velocity
        figure;
        patch(coord_x,coord_y,vel_y,'EdgeColor',edge_color);
        title('Y-Velocity');
        box on;
        axis xy equal;
        axis(axis_domain);
        set(gca,'Layer','top');
        
        % save figure
        if (save_velocityf)
            fname = [case_name,'_y-velocity_',num2str(i,'%04d')];
            print('-djpeg','-r200',[fname,'.jpg']);
            close;
        end
        
    end
    
    % just velocity grid
    if (plot_vgrid_only)
        
        figure;
        patch(reshape(grid_vnodes(1,grid_vnodes2elements(1:4,:)),4,no_el), ...
            reshape(grid_vnodes(2,grid_vnodes2elements(1:4,:)),4,no_el), ...
            0,'EdgeColor','k','FaceColor','None');
        title('Velocity Grid');
        box on;
        axis xy equal;
        axis(axis_domain);
        set(gca,'Layer','top');
        
    end
    
end

end
