%% Implementation of Matrices A, S, A_tilde and B

clear

% Grid parameters
mx = 140; % Grid points
my = 140;
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

S11 = P_W * (-(lambda_matrix + 2*mu_matrix) * Dmx_2d) + ...
      P_E * ( (lambda_matrix + 2*mu_matrix) * Dmx_2d) + ...
      P_N * ( mu_matrix * Dmy_2d) + ...
      P_S * (-mu_matrix* Dmy_2d);

S12 = (P_W * (-lambda_matrix * Dmy_2d)) + ...
      (P_E * ( lambda_matrix * Dmy_2d)) + ...
      (P_N * ( mu_matrix * Dmx_2d))     + ...
      (P_S * (-mu_matrix * Dmx_2d));

S21 = (P_W * (-mu_matrix * Dmy_2d))     + ...
      (P_E * ( mu_matrix * Dmy_2d))     + ...
      (P_N * ( lambda_matrix * Dmx_2d)) + ...
      (P_S * (-lambda_matrix * Dmx_2d));
 
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


%% Frequency domain 

freq = 9.5; % Hz
alpha = 2*pi*freq; 

% Implementation of matrix from the time independent problem formulation
M11 = -(A_tilde + Rho_Full * alpha^2); 
M12 = -alpha*B;
M21 = alpha*B;
M22 = -(A_tilde + Rho_Full * alpha^2);

M = [M11, M12;
    M21, M22];

% Point force  
A = 10e6; % Amplitude 
ix = mx * 0.05; % Coordinates for point force 
iy = my*0.95;
d_matrix = zeros(my, mx);
d_matrix(iy, ix) = 1; 
d_spatial = d_matrix(:);
d_spatial = (1./diag(H_2d)) .* d_spatial;

f_vec_j = zeros(2 * mx * my, 1);

f_vec_j(1 : mx*my) = A * d_spatial; 

f = [f_vec_j; 
    zeros(size(f_vec_j))]; 

% Solve equation system
sol = M \ f;

s1 = reshape(sol(1 : mx*my), my, mx);
s2 = reshape(sol(2 * mx*my + 1 : 3 * mx*my), my, mx);

r1 = reshape(sol(2 * mx*my + 1 : 3 * mx*my), my, mx);
r2 = reshape(sol(3 * mx*my + 1 : 4 * mx*my), my, mx);

% Amplitude 
amp_s = sqrt(s1.^2 + s2.^2);
amp_r = sqrt(r1.^2 + r2.^2);

max_amp = max(amp_s, amp_r); % Maximum amplitude

title_text = sprintf('Displacement at frequency f = %.2f Hz', freq);

file_name = sprintf('plot_no_room_mx%d_my%d_f%.2f.png', mx, my, freq);

% Plot
fig = figure;
imagesc(x, y, max_amp);
axis xy; 
axis square;
colorbar;
clim([0, max(max_amp(:)) * 1]);
xlabel('Domain width [meter]');
ylabel('Domain height [meter]');
title(title_text);
saveas(fig, file_name);

disp('Code finished')