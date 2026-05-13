%% --- Implementation without SAT ---

clear

% Grid parameters
mx = 200; % Grid points
my = 200;
Lx = 10; % Domain lenght
Ly = 10;
hx = Lx/(mx-1); % Step lenght 
hy = Ly/(my-1);
margin = round(mx*0.2);

x = linspace(0, Lx, mx).' ; % Define x and y vector 
y = linspace(0, Ly, my).' ;
[X, Y] = meshgrid(x, y); % Create grid in 2D

% Material constants
lambda = 42.3 * 10^9;
mu = 22.7 * 10^9;
rho = 2700;

% Create 1D operators for x and y 
[Dpx, Dmx, Hx] = SBP_5(mx, hx);
[Dpy, Dmy, Hy] = SBP_5(my, hy);

% Create identity matrices 
Ix = eye(mx); 
Iy = eye(my);

% Create 2D operatoors for x, y
Dpx_2d = kron(sparse(Dpx), sparse(Iy)); % Derivative in x direction with Dp
Dpy_2d = kron(sparse(Ix), sparse(Dpy)); % Derivative in y direction wuth Dp

Dmx_2d = kron(sparse(Dmx), sparse(Iy)); % Derivative in x direction with Dm
Dmy_2d = kron(sparse(Ix), sparse(Dmy)); % Derivative in x direction with Dm

% Create H for 2D
H_2d = kron(sparse(Hx), sparse(Hy));

% Elements in derivative matrix A
A11 = Dpx_2d * (lambda + 2*mu) * Dmx_2d + Dpy_2d * mu * Dmy_2d;
A12 = Dpy_2d * mu * Dmx_2d + Dpx_2d * lambda * Dmy_2d; 
A21 = Dpx_2d * mu * Dmy_2d + Dpy_2d * lambda * Dmx_2d;
A22 = Dpx_2d * mu * Dmx_2d + Dpy_2d * (lambda + 2*mu) * Dmy_2d;

A = [A11, A12;
    A21, A22];

%%

% SAT

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


% Create parts of S matrix

% Block S11: Contribution to Fx from vx
S11 = sparse((H_2d \ e_W * Hy * e_W') * (-(lambda + 2*mu) * Dmx_2d) + ...
      (H_2d \ e_E * Hy * e_E') * ( (lambda + 2*mu) * Dmx_2d) + ...
      (H_2d \ e_N * Hx * e_N') * ( mu * Dmy_2d) + ...
      (H_2d \ e_S * Hx * e_S') * (-mu * Dmy_2d));

% Block S12: Contribution to Fx from vy
S12 = sparse((H_2d \ e_W * Hy * e_W') * (-lambda * Dmy_2d) + ...
      (H_2d \ e_E * Hy * e_E') * ( lambda * Dmy_2d) + ...
      (H_2d \ e_N * Hx * e_N') * ( mu * Dmx_2d) + ...
      (H_2d \ e_S * Hx * e_S') * (-mu * Dmx_2d));

% Block S21: Contribution to Fy from vx 
S21 = sparse((H_2d \ e_W * Hy * e_W') * (-mu * Dmy_2d) + ...
      (H_2d \ e_E * Hy * e_E') * ( mu * Dmy_2d) + ...
      (H_2d \ e_N * Hx * e_N') * ( lambda * Dmx_2d) + ...
      (H_2d \ e_S * Hx * e_S') * (-lambda * Dmx_2d));

% Block S22: Contribution to Fy from vy 
S22 = sparse((H_2d \ e_W * Hy * e_W') * (-mu * Dmx_2d) + ...
      (H_2d \ e_E * Hy * e_E') * ( mu * Dmx_2d) + ...
      (H_2d \ e_N * Hx * e_N') * ( (lambda + 2*mu) * Dmy_2d) + ...
      (H_2d \ e_S * Hx * e_S') * (-(lambda + 2*mu) * Dmy_2d));

% Total SAT matrix
S = [S11, S12; 
     S21, S22];

A_tilde = A - S; 


% construct matrix B in v_tt = A_tilde*v + B*v_t

% p and s wave velocities
cp = sqrt((lambda + 2*mu) / rho);
cs = sqrt(mu / rho);

% Impedances
Zp = rho * cp;
Zs = rho * cs;

% All North-terms zero due to traction free BC:s
B11 = sparse((H_2d \ e_W * Hy * e_W') * (-Zp) + ...
      (H_2d \ e_E * Hy * e_E') * (-Zp) + ...
      (H_2d \ e_S * Hx * e_S') * (-Zs));

% No contribution from v_dot_y (y velocity) in the equation for x
B12 = sparse(size(B11,1), size(B11,2));

B22 = sparse((H_2d \ e_W * Hy * e_W') * (-Zs) + ...
      (H_2d \ e_E * Hy * e_E') * (-Zs) + ...
      (H_2d \ e_S * Hx * e_S') * (-Zp));

% No contribution from v_dot_x (x velocity) in the equation for y
B21 = sparse(size(B22,1), size(B22,2));

%total matrix B
B = [B11, B12; 
     B21, B22];


%% RK4 

% Time 
dt = 0.00001;          % Time step 
T_final = 0.002;       % Total time
t_steps = T_final/dt;

% Gaussian pulse in time
t0 = 0.0005;      
sigma_t = 0.00001; % width

% f(t) 
A = 1e6;
f_t = @(t) A * exp(-(t-t0)^2 / (2*sigma_t^2));

% point force
ix = round(mx/2); 
iy = round(my/2); 
node_idx = (ix-1)*my + iy; % location for point force in 2d flattened vector

d_spatial = zeros(mx*my, 1);
d_spatial(node_idx) = 1; 

d_spatial = (1./diag(H_2d)) .* d_spatial;

d_full = zeros(4*mx*my, 1);
d_full(node_idx) = d_spatial(node_idx);


v_initial = zeros(2*mx*my, 1);
w_initial = zeros(2*mx*my, 1);
u = [v_initial; w_initial];

% Definition of function f(u, t)
f = @(u, t) [u(2*mx*my+1:end); 
             (1/rho) * (A_tilde * u(1:2*mx*my) + B * u(2*mx*my+1:end) + d_full(1:2*mx*my) * f_t(t))];

             
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
        v_x_max = max(abs(u(2*mx*my:3*mx*my))); 
        fprintf('Time: %.5f s | Maxvelocity: %e m/s\n', t, v_x_max);
        v_current = u(1:mx*my);
        max_v = max(abs(v_current));
        V_plot = reshape(v_current, my, mx);
        imagesc(x, y, V_plot); 
        axis xy;
        colorbar;
        title(['Time: ', num2str(n*dt)]);
        drawnow;
    end
end

