function sdvig
% $Id$

clear all;
close all;
clc;

% name of the case
case_name = 'extension';

% perform output?
output = true;

% check if it's restart
last_step = csvread('nstep');
if (last_step == 0)
    restart = false;
else
    restart = true;
    restart_file = [case_name,'_',num2str(last_step,'%04d'),'.h5'];
end

% plasticity
plasticity = true;
pl_iter_max = 15;
plast_tol = 1.0e-3;
% Drucker-Prager
phi = atan(0.0);
cohesion_0 = 4.0;
cohesion_1 = 1.0;
cosphi = cos(phi);
sinphi = sin(phi);
strain_min = 0.0;
strain_max = 0.05;
strain_dif = 1.0 / (strain_max - strain_min);

% mesh and domain
no_el_x = 192;
no_el_y = 64;
xmin0 = -1.5;
xmax0 =  1.5;
ymin0 =  0.0;
ymax0 =  1.0;

% viscosity
nu = 100.0;
inclusion_nu = 1.0;
min_nu = 0.01;

% density
rho = 1.0;

% air layer
airl = 12;
air_nu = min_nu;
air_rho = 0.0;

% time parameters
ttotal  = 40.0;
dt      =  0.1;
ittotal = ceil(ttotal / dt);
if (restart)
    itstart = hdf5read(restart_file,'/Info/CurrentIteration') + 1;
else
    itstart = 1;
end

% external force field
G = [0.0 -10.0];

% boundary speed
side_u = 0.05;

% deform mesh (according to boundary speed) ?
mesh_def = true;

% deformed mesh
if (restart)
    domain = hdf5read(restart_file,'/Info/Domain');
    xmin = domain(1);
    xmax = domain(2);
    ymin = domain(3);
    ymax = domain(4);
else
    xmin = xmin0;
    xmax = xmax0;
    ymin = ymin0;
    ymax = ymax0;    
end

% form of D matrix (if true - use compressible, incompressible otherwise)
D_comp = false;

% D matrix
if (D_comp)
    D = [ 4/3 -2/3 0; ...
         -2/3  4/3 0; ...
            0    0 1];
else
    D = [ 2  0  0; ...
          0  2  0; ...
          0  0  1];
end

% element type (1 - Q1P0, 2 - Q1Q1, 3 - Q2P-1)
el_type = 3;

switch el_type
    case 1
        % Q1P0 linearV-constantP(discontinuously)
        no_vnodes_el = 4;
        no_pnodes_el = 1;
        no_vnodes_x = no_el_x + 1;
        no_vnodes_y = no_el_y + 1;
        no_pnodes_x = no_el_x;
        no_pnodes_y = no_el_y;
        vxvec = linspace(xmin,xmax,no_vnodes_x);
        vyvec = linspace(ymin,ymax,no_vnodes_y);
        pxvec = 0.5*(vxvec(1:end-1) + vxvec(2:end));
        pyvec = 0.5*(vyvec(1:end-1) + vyvec(2:end));
        [vX, vg_num, vbpoints] = generate_mesh({vxvec vyvec}, ...
            [no_el_x no_el_y], [no_vnodes_x no_vnodes_y], no_vnodes_el);
        [pX, pg_num, pbpoints] = generate_mesh({pxvec pyvec}, ...
            [no_el_x no_el_y], [no_pnodes_x no_pnodes_y], no_pnodes_el);
        no_el = no_el_x * no_el_y;
        no_vnodes = no_vnodes_x * no_vnodes_y;
        no_pnodes = no_pnodes_x * no_pnodes_y;
        nip = 4;
    case 2
        % Q1Q1 linearV-linearP(continuously)
        no_vnodes_el = 4;
        no_pnodes_el = 4;
        no_vnodes_x = no_el_x + 1;
        no_vnodes_y = no_el_y + 1;
        no_pnodes_x = no_vnodes_x;
        no_pnodes_y = no_vnodes_y;
        vxvec = linspace(xmin,xmax,no_vnodes_x);
        vyvec = linspace(ymin,ymax,no_vnodes_y);
        pxvec = vxvec;
        pyvec = vyvec;
        [vX, vg_num, vbpoints] = generate_mesh({vxvec vyvec}, ...
            [no_el_x no_el_y], [no_vnodes_x no_vnodes_y], no_vnodes_el);
        [pX, pg_num, pbpoints] = generate_mesh({pxvec pyvec}, ...
            [no_el_x no_el_y], [no_pnodes_x no_pnodes_y], no_pnodes_el);
        no_el = no_el_x * no_el_y;
        no_vnodes = no_vnodes_x * no_vnodes_y;
        no_pnodes = no_pnodes_x * no_pnodes_y;
        nip = 4;
        alpha = 1.0;
    case 3
        % Q2P-1 quadraticV-linearP(discontinuously)
        no_vnodes_el = 9;
        no_pnodes_el = 1;
        no_vnodes_x = no_el_x * 2 + 1;
        no_vnodes_y = no_el_y * 2 + 1;
        no_pnodes_x = no_el_x;
        no_pnodes_y = no_el_y;
        vxvec = linspace(xmin,xmax,no_vnodes_x);
        vyvec = linspace(ymin,ymax,no_vnodes_y);
        pxvec = vxvec(2:2:end-1);
        pyvec = vyvec(2:2:end-1);
        [vX, vg_num, vbpoints] = generate_mesh({vxvec vyvec}, ...
            [no_el_x no_el_y], [no_vnodes_x no_vnodes_y], no_vnodes_el);
        [pX, pg_num, pbpoints] = generate_mesh({pxvec pyvec}, ...
            [no_el_x no_el_y], [no_pnodes_x no_pnodes_y], no_pnodes_el);
        no_el = no_el_x * no_el_y;
        no_vnodes = no_vnodes_x * no_vnodes_y;
        no_pnodes = no_pnodes_x * no_pnodes_y;
        % add nodes for pressure gradients (TODO: move it to generate_mesh)
        clear pg_num;
        pg_num = reshape([1:3*no_pnodes],3,no_pnodes); % !!! pX is not consistent anymore
        pg_num = uint32(pg_num);
        no_pnodes = no_pnodes * 3;
        no_pnodes_el = 3;
        nip = 9;
end

% viscosity and density
ERho = zeros(no_el,1);
ENu = zeros(no_el,1);

% coordinates of elements (for particles search)
el_xvec = linspace(xmin,xmax,no_el_x+1);
el_yvec = linspace(ymin,ymax,no_el_y+1);
% particles' data
if (restart)
    % read data
    part_per_side = hdf5read(restart_file,'/Info/ParticlesPerSide');
    XPart = hdf5read(restart_file,'/Particles/XCoordinates');
    YPart = hdf5read(restart_file,'/Particles/YCoordinates');
    no_part = size(XPart,1);
    RhoPart = hdf5read(restart_file,'/Particles/Density');
    NuPart = hdf5read(restart_file,'/Particles/Viscosity');
    ColorPart = hdf5read(restart_file,'/Particles/Material');
    StrainPl = hdf5read(restart_file,'/Particles/PlasticStrain');
else
    % exclude air layer
    no_el_yp = no_el_y - airl;
    ymaxp = el_yvec(end-airl);
    % create particles
    part_per_side = 3;
    no_part_x = no_el_x * part_per_side;
    no_part_y = no_el_yp * part_per_side;
    no_part = no_part_x * no_part_y;
    dx = (xmax - xmin) / no_part_x; dx2 = dx / 2;
    dy = (ymaxp - ymin) / no_part_y; dy2 = dy / 2;
    xvec_part = linspace(xmin+dx2,xmax-dx2,no_part_x);
    yvec_part = linspace(ymin+dy2,ymaxp-dy2,no_part_y);
    XPart = reshape(repmat(xvec_part,no_part_y,1),1,no_part)';
    YPart = repmat(yvec_part,1,no_part_x)';
    c = 0.0; % insert some noise
    rand('twister',0); % init number generator
    XPart = XPart + c * dx * (rand(no_part,1) - 0.5);
    YPart = YPart + c * dy * (rand(no_part,1) - 0.5);
    % viscosity and density
    RhoPart = zeros(no_part,1);
    NuPart = zeros(no_part,1);
    RhoPart(:) = rho;
    NuPart(:) = nu;
    % insert inclusion
    for ip = 1:no_part
        x = XPart(ip);
        y = YPart(ip);
        if ((x > -0.02 && x < 0.02) && (y < 0.04))
            NuPart(ip)  = inclusion_nu;
        end
    end
    % paint particles (for visualization)
    ColorPart = zeros(no_part,1);
    lthick = 4; % layer thickness (number of elements)
    color1 = 2;
    color2 = 3;
    for ip = 1:no_part
        y = YPart(ip);
        no_layer = fix(y / el_yvec(1+lthick));
        if (~rem(no_layer,2))
            % even layer
            ColorPart(ip) = color1;
        else
            % odd layer
            ColorPart(ip) = color2;
        end
    end
    % save space
    ColorPart = uint8(ColorPart);
    % plastic strain
    StrainPl = zeros(no_part,1);
end
% velocities
VxPart = zeros(no_part,1);
VyPart = zeros(no_part,1);
% 2nd invariant strain rate
Strain2 = zeros(no_part,1);

% BC
bc = true;
% velocity boundary points
vleft    = vbpoints{1};
vright   = vbpoints{2};
vtop     = vbpoints{3};
vbottom  = vbpoints{4};
lvleft   = length(vleft);
lvright  = length(vright);
lvtop    = length(vtop);
lvbottom = length(vbottom);
% x-velocity
vbcdofx = [vleft vright]';
vbcvalx = zeros(length(vbcdofx),1);
vbcvalx(1:lvleft)     = -side_u; % left
vbcvalx(lvleft+1:end) =  side_u; % right
% y-velocity
vbcdofy = [vbottom]';
vbcvaly = zeros(length(vbcdofy),1);
vbcvaly(1:lvbottom) = 0.0; % bottom

% DOF velocity
vdof_node = 2;                    % per node
vdof_el = vdof_node*no_vnodes_el; % per element
vdof = vdof_node*no_vnodes;       % total

% DOF pressure
pdof_node = 1;                    % per node
pdof_el = pdof_node*no_pnodes_el; % per element
pdof = pdof_node*no_pnodes;       % total

% % DOF total
% dof_el = vdof_el + pdof_el;
% dof = vdof + pdof;

% aux. arrays
vn_f = reshape((1:vdof),vdof_node,no_vnodes);     % node -> eq.numbers
pn_f = reshape((1:pdof),pdof_node,no_pnodes);     % node -> eq.numbers
vg_g = zeros(vdof_el,no_el);                      % element -> eq.numbers
pg_g = zeros(pdof_el,no_el);                      % element -> eq.numbers
for iel = 1:no_el
    neq = vn_f(:,vg_num(:,iel));
    vg_g(:,iel) = reshape(neq,vdof_el,1);
    neq = pn_f(:,pg_num(:,iel));
    pg_g(:,iel) = reshape(neq,pdof_el,1);
end

% integration points and weights
[ipx, ipw] = ips(nip);

% shape functions and their derivatives wrt local coordinates
[Nv, dNvdu] = shp_deriv(ipx, no_vnodes_el, true); % velocity
Np = shp_deriv(ipx, no_pnodes_el, false); % pressure

% allocate matrices and vectors
N_A_triplets = no_el*vdof_el*vdof_el;
I_A = zeros(N_A_triplets,1);
J_A = zeros(N_A_triplets,1);
X_A = zeros(N_A_triplets,1);
N_Q_triplets = no_el*vdof_el*pdof_el;
I_Q = zeros(N_Q_triplets,1);
J_Q = zeros(N_Q_triplets,1);
X_Q = zeros(N_Q_triplets,1);
N_M_triplets = no_el*pdof_el*pdof_el;
I_M = zeros(N_M_triplets,1);
J_M = zeros(N_M_triplets,1);
X_M = zeros(N_M_triplets,1);
X_C = zeros(N_M_triplets,1);
RHS = zeros(vdof,1);
A_loc = zeros(vdof_el,vdof_el);
A_loc_all = cell(no_el,1);
Q_loc = zeros(vdof_el,pdof_el);
M_loc = zeros(pdof_el,pdof_el);
C_loc = zeros(pdof_el,pdof_el);
QMQ_loc_all = cell(no_el,1);
RHS_loc = zeros(vdof_el,1);
B = zeros(3,vdof_el);
u = zeros(vdof,1);
uprev = zeros(vdof,1);
p = zeros(pdof,1);

% iterate over time
for itcur = itstart:1:ittotal
        
    % current time
    tcur = itcur * dt;
    
    % be verbose
    msg = sprintf('iter: %5d        time: %10.5f', itcur, tcur);
    disp(msg);
    
    % assign particles to elements
    clear Part2El;
    Part2El = cell(no_el,1);
    for ip = 1:no_part
        i = find(el_xvec < XPart(ip),1,'last');
        j = find(el_yvec < YPart(ip),1,'last');
        % check if particle is inside the domain
        if (~isempty(i) && i < (no_el_x+1) && ...
                ~isempty(j) && j < (no_el_y+1))
            iel = (i-1)*no_el_y + j;
            Part2El{iel} = [Part2El{iel} ip];
        end
    end
    
    % particles->mesh interpolation
    for iel = 1:no_el
        % particles inside the element
        epart = Part2El{iel};
        % empty element
        if (isempty(epart))
            ENu(iel) = air_nu;
            ERho(iel) = air_rho;
        else
            % arithmetic averaging
            ENu(iel) = mean(NuPart(epart));
            ERho(iel) = mean(RhoPart(epart));
        end
    end
    
    % init matrices
    RHS(:) = 0;
    A_triplets = 0;
    Q_triplets = 0;
    M_triplets = 0;
    
    % incompressibility factor
    % k = 1.0e+3 * max(ENu);
    k = 1.0e+3 * nu;
    
    % build matrices
    for iel = 1:no_el
        % init local matrices
        A_loc(:) = 0;
        Q_loc(:) = 0;
        M_loc(:) = 0;
        C_loc(:) = 0;
        RHS_loc(:) = 0;
        % global coordinates
        Xi = vX(:,vg_num(:,iel))';
        % loop over integration points
        for ip = 1:nip
            Nvi = Nv(ip,:);
            Npi = Np(ip,:);
            dNvdui = dNvdu(2*ip-1:2*ip,:);
            J = dNvdui * Xi;
            detJ = det(J);
            invJ = inv(J);
            dNvdxi = invJ * dNvdui;   % derivatives wrt global coordinates
            B(1,1:2:end-1) = dNvdxi(1,:);
            B(2,2:2:end)   = dNvdxi(2,:);
            B(3,1:2:end-1) = dNvdxi(2,:);
            B(3,2:2:end)   = dNvdxi(1,:);
            Bvol = dNvdxi;
            w = ipw(ip) * detJ;
            A_loc = A_loc + w * (B' * D * B);
            Q_loc = Q_loc - w * (Bvol(:) * Npi);
            if (el_type == 1 || el_type == 3)
                % mass matrix for Q1P0 or Q2P-1 elements
                M_loc = M_loc + w * (Npi' * Npi);
            else
                % stabilization matrix for Q1Q1
                C_loc = C_loc + w * ((Npi - 0.25)' * (Npi - 0.25));
            end
            r = G' * Nvi;
            RHS_loc = RHS_loc + w * r(:);
        end
        
        % store A_loc
        A_loc_all{iel} = A_loc;
        
        if (el_type == 1 || el_type == 3)
            % Q1P0 or Q2P-1 elements (discontinuous pressure)
                        
            % static condensation (TODO: check)
            invM_loc = inv(M_loc);
            QMQ_loc = Q_loc * invM_loc * Q_loc';
            A_loc = ENu(iel) * A_loc + k * QMQ_loc;            
            
            % store QMQ_loc
            QMQ_loc_all{iel} = QMQ_loc;
        else
            % Q1Q1 element
            
            A_loc = ENu(iel) * A_loc;
            C_loc = (alpha / ENu(iel)) * C_loc;
        end
        
        % construct triplets
        vn = vg_g(:,iel);
        pn = pg_g(:,iel);
        for i = 1:vdof_el
            for j = 1:vdof_el
                A_triplets = A_triplets + 1;
                I_A(A_triplets) = vn(i);
                J_A(A_triplets) = vn(j);
                X_A(A_triplets) = A_loc(i,j);
            end
            for j = 1:pdof_el
                Q_triplets = Q_triplets + 1;
                I_Q(Q_triplets) = vn(i);
                J_Q(Q_triplets) = pn(j);
                X_Q(Q_triplets) = Q_loc(i,j);
            end
        end
        for i = 1:pdof_el
            for j = 1:pdof_el
                M_triplets = M_triplets + 1;
                I_M(M_triplets) = pn(i);
                J_M(M_triplets) = pn(j);
                if (el_type == 1 || el_type == 3)
                    % Q1P0 or Q2P-1
                    X_M(M_triplets) = invM_loc(i,j);
                else
                    % Q1Q1
                    X_C(M_triplets) = C_loc(i,j);
                end
            end
        end
        
        % update RHS
        RHS(vn) = RHS(vn) + ERho(iel) * RHS_loc;
        
    end

    % construct matrices from triplets
    Q = sparse(I_Q,J_Q,X_Q,vdof,pdof);
    if (el_type == 1 || el_type == 3)
        % Q1P0 or Q2P-1
        invM = sparse(I_M,J_M,X_M,pdof,pdof);
    else
        % Q1Q1
        C = sparse(I_M,J_M,X_C,pdof,pdof);
    end                                        
    
    % for plastic iterations
    uprev(:) = 0;
    pl_iter = 0;
    
    while (true)
        
        % matrix A
        A = sparse(I_A,J_A,X_A,vdof,vdof);
        
        % impose BC
        if (bc)
            % equations corresponding to BC
            neqx = vn_f(1,vbcdofx);
            neqy = vn_f(2,vbcdofy);
            neq = [neqx neqy];
            % BC values
            bcval = [vbcvalx; vbcvaly];
            u(neq) = bcval;
            % free equations
            free_eq = [1:1:vdof];
            free_eq(neq) = [];
            if (el_type == 1 || el_type == 3)
                % Q1P0 or Q2P-1
                % impose BC on RHS
                RHSs = RHS - A(:,neq) * bcval;
                % exclude BC equations from A
                A = A(free_eq,free_eq);
            else
                % Q1Q1
                % impose BC on RHS for momentum eq.
                RHSs = RHS - A(:,neq) * bcval;
                % impose BC on RHS for continuity eq.
                Qt = Q';
                RHSs1 = zeros(size(p)) - Qt(:,neq) * bcval;
                % combine
                RHSs = [RHSs(free_eq);RHSs1];
                % exclude BC equations from A and Q
                A = A(free_eq,free_eq);
                Qf = Q(free_eq,:);
            end
        end
        
        if (el_type == 1 || el_type == 3)
            % Q1P0 or Q2P-1
            % PH iterations
            p(:) = 0;
            j = 0;
            divmax = 1;
            while ((divmax > 1e-10) && (j < 10))
                j = j + 1;
                RHSss = RHSs - Q * p;
                % velocity
                u(free_eq) = A \ RHSss(free_eq);
                % divergence
                div = invM * (Q' * u); % TODO: can it be done like static condensation (invM*Q') ?
                divmax = max(abs(div(:)));
                % pressure
                p = p + k * div;
            end
        else
            % Q1Q1
            A = [A, Qf; Qf', -C];
            x = A \ RHSs;
            n = length(free_eq);
            u(free_eq) = x(1:n);
            p = x(n+1:end);
        end
        
        % plasticity
        if (plasticity)
            
            % residue
            pl_res = max(abs(u - uprev)) / max(abs(u));
            
            % be verbose
            msg = sprintf('plasticity: %5d  %20.10f', pl_iter, pl_res);
            disp(msg);
            
            % check convergence
            if (pl_res < plast_tol || pl_iter == pl_iter_max)
                break;
            end
            
            % iteration number
            pl_iter = pl_iter + 1;
            
            % matrix A
            clear A;
            A_triplets = 0;
            
            % loop over elements, correct viscosity and update matrix A
            for iel = 1:no_el
                % particles inside the element
                epart = Part2El{iel};
                % skip empty elements
                if (~isempty(epart))
                    % global coordinates
                    Xi = vX(:,vg_num(:,iel));
                    elX = Xi(1,:);
                    elY = Xi(2,:);
                    % velocities
                    uel = u(vg_g(:,iel));
                    % pressure
                    pel = p(pg_g(:,iel));
                    % loop over particles inside the element
                    for ip = epart
                        % local particle coordinates
                        XiP = ((XPart(ip) - elX(1)) / (elX(2) - elX(1))) * 2 - 1;
                        EtaP = ((YPart(ip) - elY(2)) / (elY(3) - elY(2))) * 2 - 1;
                        % compute strain rate at particle position
                        dNdui = deriv_p([XiP EtaP], no_vnodes_el);
                        J = dNdui * Xi';
                        dNdxi = inv(J) * dNdui; % wrt global coordinates
                        B(1,1:2:end-1) = dNdxi(1,:);
                        B(2,2:2:end)   = dNdxi(2,:);
                        B(3,1:2:end-1) = dNdxi(2,:);
                        B(3,2:2:end)   = dNdxi(1,:);
                        Strain = B * uel;
                        Strain(3) = 0.5 * Strain(3);
                        % 2nd invariant strain rate
                        Strain2(ip) = sqrt(0.5 * (Strain(1)^2 + ...
                            Strain(2)^2) + Strain(3)^2);
                        % 2nd invariant stress
                        Stress2 = 2.0 * nu * Strain2(ip);                        
                        % pressure
                        switch el_type
                            case 1
                                % Q1P0
                                pp = pel;
                            case 2
                                % Q1Q1
                                Ni = shp_p([XiP EtaP], no_pnodes_el);
                                pp = Ni * pel;
                            case 3
                                % Q2P-1
                                pp = pel(1) + XiP * pel(2) + EtaP * pel(3);
                        end
                        % strain weakening
                        if (StrainPl(ip) > strain_max)
                            % maximum weakening
                            cohesion = cohesion_1;
                        elseif (StrainPl(ip) < strain_min)
                            % no weakening
                            cohesion = cohesion_0;
                        else
                            % linear weakening with strain
                            cohesion = cohesion_0 + (cohesion_1 - cohesion_0) * ...
                                (StrainPl(ip) - strain_min) * strain_dif;
                        end
                        % yielding stress
                        Syield = cosphi * cohesion + sinphi * pp;
                        % check yielding
                        if (Stress2 - Syield > 0)
                            % compute effective viscosity
                            NuPart(ip) = 0.5 * Syield / Strain2(ip);
                        else
                            % set bgrnd viscosity (TODO: does it needed?)
                            NuPart(ip) = nu;
                        end
                        % cutoff viscosity
                        if (NuPart(ip) < min_nu)
                            NuPart(ip) = min_nu;
                        end
                    end
                    % arithmetic averaging
                    ENu(iel) = mean(NuPart(epart));
                end
                
                % update matrix A
                vn = vg_g(:,iel);
                if (el_type == 1 || el_type == 3)
                    % Q1P0 or Q2P-1
                    A_loc = ENu(iel) * A_loc_all{iel} + k * QMQ_loc_all{iel};
                else
                    % Q1Q1
                    A_loc = ENu(iel) * A_loc_all{iel};
                end
                for i = 1:vdof_el
                    for j = 1:vdof_el
                        A_triplets = A_triplets + 1;
                        I_A(A_triplets) = vn(i);
                        J_A(A_triplets) = vn(j);
                        X_A(A_triplets) = A_loc(i,j);
                    end
                end
                
            end
            
            % save velocity
            uprev = u;
            
        else
            break;
        end
        
    end
    
    % clear matrices
    clear A;
    clear Q;
    clear invM;
    clear C
    
    % mesh->particles interpolation
    ux = u(1:2:end-1);
    uy = u(2:2:end);
    for iel = 1:no_el
        % particles inside the element
        epart = Part2El{iel};
        % empty element
        if (isempty(epart))
            continue;
        end
        % element nodes
        en = vg_num(:,iel);
        % global coordinates
        Xi = vX(:,en);
        elX = Xi(1,:);
        elY = Xi(2,:);
        % loop over particles inside the element
        for ip = epart
            % local particle coordinates
            XiP = ((XPart(ip) - elX(1)) / (elX(2) - elX(1))) * 2 - 1;
            EtaP = ((YPart(ip) - elY(2)) / (elY(3) - elY(2))) * 2 - 1;
            % interpolate velocities from nodes to particle
            Ni = shp_p([XiP EtaP], no_vnodes_el);
            VxPart(ip) = Ni * ux(en);
            VyPart(ip) = Ni * uy(en);
            % update plastic strain
            if (plasticity)
                StrainPl(ip) = StrainPl(ip) + ...
                    dt * Strain2(ip) * (1.0 - NuPart(ip) / nu);
            end
        end
    end
    
    % advect particles
    XPart = XPart + dt * VxPart;
    YPart = YPart + dt * VyPart;
    
    % deform mesh
    if (mesh_def)
        xmin = xmin - dt * side_u;
        xmax = xmax + dt * side_u;
        vxvec = linspace(xmin,xmax,no_vnodes_x);
        [vX, vg_num, vbpoints] = generate_mesh({vxvec vyvec}, ...
            [no_el_x no_el_y], [no_vnodes_x no_vnodes_y], no_vnodes_el);
        el_xvec = linspace(xmin,xmax,no_el_x+1);
    end
    
    % output
    if (output)
        fname = [case_name,'_',num2str(itcur,'%04d'),'.h5'];
        hdf5write(fname,'/Info/CurrentIteration',itcur);
        hdf5write(fname,'/Info/CurrentTime',tcur,'writemode','append');
        hdf5write(fname,'/Info/InitialDomain',[xmin0 xmax0 ymin0 ymax0],'writemode','append');
        hdf5write(fname,'/Info/Domain',[xmin xmax ymin ymax],'writemode','append');
        hdf5write(fname,'/Info/InitialResolution',[no_el_x no_el_y],'writemode','append');
        hdf5write(fname,'/Info/ElementType',el_type,'writemode','append');        
        hdf5write(fname,'/Info/ParticlesPerSide',part_per_side,'writemode','append');
        hdf5write(fname,'/Grid/VelocityNodes',vX,'writemode','append');
        hdf5write(fname,'/Grid/VelocityNodes2Elements',vg_num,'writemode','append');
        hdf5write(fname,'/Grid/XVelocity',ux,'writemode','append');
        hdf5write(fname,'/Grid/YVelocity',uy,'writemode','append');
        hdf5write(fname,'/Grid/PressureNodes',pX,'writemode','append');
        hdf5write(fname,'/Grid/PressureNodes2Elements',pg_num,'writemode','append');
        hdf5write(fname,'/Grid/Pressure',p,'writemode','append');
        hdf5write(fname,'/Particles/XCoordinates',XPart,'writemode','append');
        hdf5write(fname,'/Particles/YCoordinates',YPart,'writemode','append');
        hdf5write(fname,'/Particles/XVelocity',VxPart,'writemode','append');
        hdf5write(fname,'/Particles/YVelocity',VyPart,'writemode','append');
        hdf5write(fname,'/Particles/StrainInvariant2',Strain2,'writemode','append');
        hdf5write(fname,'/Particles/PlasticStrain',StrainPl,'writemode','append');
        hdf5write(fname,'/Particles/Viscosity',NuPart,'writemode','append');
        hdf5write(fname,'/Particles/Density',RhoPart,'writemode','append');
        hdf5write(fname,'/Particles/Material',ColorPart,'writemode','append');
    end
    
    % save iteration number
    csvwrite('nstep',itcur);
    
end

end
