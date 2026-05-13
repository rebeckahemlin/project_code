%% --- Implementation without SAT ---

clear

% Grid parameters
mx = 120; % Grid points
my = 120;
Lx = 1000; % Domain lenght
Ly = 1000;
hx = Lx/(mx-1); % Step lenght 
hy = Ly/(my-1);
margin = round(mx*0.2);

x = linspace(0, Lx, mx).' ; % Define x and y vector 
y = linspace(0, Ly, my).' ;
[X, Y] = meshgrid(x, y); % Create grid i 2D


% Material values
% Parameter varlues for granite
lambda = ones(my, mx) * (42.3 * 10^9); 
mu = ones(my, mx) * (22.7 * 10^9); 
rho = ones(my, mx) * (2700); 

room_x = round(mx*0.80):round(mx*0.83);
room_y = round(my*0.977):round(my*0.997);

lambda(room_y, room_x) = 2.2 * 10^9; % Bulk modulus for water (approx. 2.2 GPa)
mu(room_y, room_x) = 0; % Shear modulus of a fluid is zero (no S-waves)
rho(room_y, room_x) = 1000; % Density of water 

% Flatten the 2D material matrices into 1D column vectors 
lambda_vec = lambda(:); 
mu_vec = mu(:); 
rho_vec = rho(:);

% Put the flattened vectors onto the diagonal
lambda_matrix = spdiags(lambda_vec, 0, mx*my, mx*my); 
mu_matrix= spdiags(mu_vec, 0, mx*my, mx*my);
rho_matrix = spdiags(rho_vec, 0, mx*my, mx*my); 
inv_rho_matrix = spdiags(1./rho_vec, 0, mx*my, mx*my);
Inv_Rho_Full = blkdiag(inv_rho_matrix, inv_rho_matrix); % for u = [vx, vy], size: (2*mx*my) x (2*mx*my)
Rho_Full = blkdiag(rho_matrix, rho_matrix);

% 1D operators for x and y 
[Dpx, Dmx, Hx] = SBP_9(mx, hx);
[Dpy, Dmy, Hy] = SBP_9(my, hy);

% identity matrices 
Ix = eye(mx); 
Iy = eye(my);

% 2D operatoors for x, y
Dpx_2d = kron(sparse(Dpx), sparse(Iy)); % Derivative in x direction with Dp
Dpy_2d = kron(sparse(Ix), sparse(Dpy)); % Derivative in y direction wuth Dp

Dmx_2d = kron(sparse(Dmx), sparse(Iy)); % Derivative in x direction with Dm
Dmy_2d = kron(sparse(Ix), sparse(Dmy)); % Derivative in x direction with Dm

% H for 2D
H_2d = kron(sparse(Hx), sparse(Hy));

% Elements in derivative matrix A
A11 = Dpx_2d * (lambda_matrix + 2*mu_matrix) * Dmx_2d + Dpy_2d * mu_matrix * Dmy_2d;
A12 = Dpy_2d * mu_matrix * Dmx_2d + Dpx_2d * lambda_matrix * Dmy_2d; 
A21 = Dpx_2d * mu_matrix * Dmy_2d + Dpy_2d * lambda_matrix * Dmx_2d;
A22 = Dpx_2d * mu_matrix * Dmx_2d + Dpy_2d * (lambda_matrix + 2*mu_matrix) * Dmy_2d;

A = [A11, A12;
    A21, A22];

%% SAT

% Create 1D boundary vectors in x and y direction 
e1x = sparse(mx, 1); e1x(1) = 1;  % Size (mx x 1)
emx = sparse(mx, 1); emx(mx) = 1; % Size (mx x 1)

e1y = sparse(my, 1); e1y(1) = 1;  % Size (my x 1)
emy = sparse(my, 1); emy(my) = 1; % Size (my x 1)


% Create in 2D
e_W = kron(e1x, sparse(eye(my))); % Left wall (West), size(mx*my) x my
e_E = kron(emx, sparse(eye(my))); % Right wall (East), size (mx*my) x my
e_S = kron(sparse(eye(mx)), e1y); % Bottom wall (South), size (mx*my) x mx
e_N = kron(sparse(eye(mx)), emy); % Top wall (North), size (mx*my) x mx

Hinv = 1 ./ full(diag(H_2d));
       
% These are mx*my x mx*my but only have non-zeros on the diagonal
Boundary_W = sparse(1:mx*my, 1:mx*my, e_W * (diag(Hy) .* ones(my, 1)), mx*my, mx*my);
Boundary_E = sparse(1:mx*my, 1:mx*my, e_E * (diag(Hy) .* ones(my, 1)), mx*my, mx*my);
Boundary_S = sparse(1:mx*my, 1:mx*my, e_S * (diag(Hx) .* ones(mx, 1)), mx*my, mx*my);
Boundary_N = sparse(1:mx*my, 1:mx*my, e_N * (diag(Hx) .* ones(mx, 1)), mx*my, mx*my);

P_W = Hinv .* Boundary_W;
P_E = Hinv .* Boundary_E;
P_S = Hinv .* Boundary_S;
P_N = Hinv .* Boundary_N;

% Block S11: Contribution to Fx from vx 
S11 = P_W * (-(lambda_matrix + 2*mu_matrix) * Dmx_2d) + ...
      P_E * ( (lambda_matrix + 2*mu_matrix) * Dmx_2d) + ...
      P_N * ( mu_matrix * Dmy_2d) + ...
      P_S * (-mu_matrix* Dmy_2d);

% Block S12: Contribution to Fx from vy 
S12 = (P_W * (-lambda_matrix * Dmy_2d)) + ...
      (P_E * ( lambda_matrix * Dmy_2d)) + ...
      (P_N * ( mu_matrix * Dmx_2d))     + ...
      (P_S * (-mu_matrix * Dmx_2d));

% Block S21: Contribution to Fy from vx 
S21 = (P_W * (-mu_matrix * Dmy_2d))     + ...
      (P_E * ( mu_matrix * Dmy_2d))     + ...
      (P_N * ( lambda_matrix * Dmx_2d)) + ...
      (P_S * (-lambda_matrix * Dmx_2d));

% Block S22: Contribution to Fy from vy 
S22 = (P_W * (-mu_matrix * Dmx_2d))             + ...
      (P_E * ( mu_matrix * Dmx_2d))             + ...
      (P_N * ( (lambda_matrix + 2*mu_matrix) * Dmy_2d)) + ...
      (P_S * (-(lambda_matrix + 2*mu_matrix) * Dmy_2d));

% Combine for A_tilde
S = [S11, S12; S21, S22];
A_tilde = A - S;


% construct matrix B in v_tt = A_tilde*v + B*v_t

% p and s wave velocities
cp = sqrt((lambda_matrix + 2*mu_matrix) / rho_matrix);
cs = sqrt(mu_matrix / rho_matrix);

cp_max = max(max(cp));

% Impedances
Zp = rho_matrix * cp;
Zs = rho_matrix * cs;

% We apply the impedance (Zp, Zs) only at the boundaries (W, E, S)
B11 = (P_W * (-Zp)) + (P_E * (-Zp)) + (P_S * (-Zs));
B22 = (P_W * (-Zs)) + (P_E * (-Zs)) + (P_S * (-Zp));

% B12 and B21 are zero matrices of the same size
B12 = sparse(mx*my, mx*my);
B21 = sparse(mx*my, mx*my);

B = [B11, B12; B21, B22];


%% RK4 

clear A S S11 S12 S21 S22 Boundary_N Boundary_S Boundary_E Boundary_W P_N P_S P_E P_W

% Time 
T_final = 10e-1;       % Total time
dt = 0.5 * min(hx, hy) / cp_max;          % Time step 
t_steps = T_final/dt;

% Gaussian pulse in time
t0 = 0.001;      
sigma_t = 0.0001; % width

% f(t) 
A = 1e9;
f_t = @(t) A * exp(-(t-t0)^2 / (2*sigma_t^2));

% point force
ix = 1; 
iy = my; % Upper left corner 

d_matrix = zeros(my, mx); 
% n = round(0.001*mx); % n times n "block" of point force
% d_matrix(iy-n:iy+n, ix-n:ix+n) = 1; 
d_matrix(iy, ix) = 1;

d_spatial = d_matrix(:); 
d_spatial = (1./diag(H_2d)) .* d_spatial; 

d_full = zeros(4*mx*my, 1);
d_full(1:mx*my) = d_spatial;


v_initial = zeros(2*mx*my, 1);
w_initial = zeros(2*mx*my, 1);
u = [v_initial; w_initial];

% Definition of function f(u, t) aka RHS
f = @(u, t) [u(2*mx*my+1:end); 
             Inv_Rho_Full * (A_tilde * u(1:2*mx*my) +  B * u(2*mx*my+1:end) + d_full(1:2*mx*my) * f_t(t))];

             
% --- RK4 Loop ---
for n = 1:t_steps
    t = (n-1)*dt;
    
    k1 = f(u, t);
    k2 = f(u + dt/2 * k1, t + dt/2);
    k3 = f(u + dt/2 * k2, t + dt/2);
    k4 = f(u + dt * k3, t + dt);
    
    u = u + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
    
    % Plot
    if mod(n, 20) == 0
        v_current = u(1:mx*my);
        max_v = max(abs(v_current));
        V_plot = reshape(v_current, my, mx);
        imagesc(x, y, V_plot);
        axis xy;
        %clim([-1e-4, 7e-4]);
        colorbar;
        title(['Time: ', num2str(n*dt)]);
        drawnow;
    end
end
