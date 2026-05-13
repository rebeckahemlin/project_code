
%% --- Implementation without SAT ---

clear

% Grid parameters
mx = 6; % Grid points
my = 10;
Lx = 10; % Domain lenght
Ly = 10;
hx = Lx/(mx-1); % Step lenght 
hy = Ly/(my-1);
margin = round(mx*0.2);

x = linspace(0, Lx, mx).' ; % Define x and y vector 
y = linspace(0, Ly, my).' ;
[X, Y] = meshgrid(x, y); % Create grid i 2D

% Material constants
lambda = 42.3 * 10^9;
mu = 22.7 * 10^9;
roh = 2700;

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

% SAT implementation

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


% SAT terms
S11 = (H_2d \ e_W * Hy * e_W') * (-(lambda + 2*mu) * Dmx_2d) + ...
      (H_2d \ e_E * Hy * e_E') * ( (lambda + 2*mu) * Dmx_2d) + ...
      (H_2d \ e_N * Hx * e_N') * ( mu * Dmy_2d) + ...
      (H_2d \ e_S * Hx * e_S') * (-mu * Dmy_2d);

S12 = (H_2d \ e_W * Hy * e_W') * (-lambda * Dmy_2d) + ...
      (H_2d \ e_E * Hy * e_E') * ( lambda * Dmy_2d) + ...
      (H_2d \ e_N * Hx * e_N') * ( mu * Dmx_2d) + ...
      (H_2d \ e_S * Hx * e_S') * (-mu * Dmx_2d);

S21 = (H_2d \ e_W * Hy * e_W') * (-mu * Dmy_2d) + ...
      (H_2d \ e_E * Hy * e_E') * ( mu * Dmy_2d) + ...
      (H_2d \ e_N * Hx * e_N') * ( lambda * Dmx_2d) + ...
      (H_2d \ e_S * Hx * e_S') * (-lambda * Dmx_2d);

S22 = (H_2d \ e_W * Hy * e_W') * (-mu * Dmx_2d) + ...
      (H_2d \ e_E * Hy * e_E') * ( mu * Dmx_2d) + ...
      (H_2d \ e_N * Hx * e_N') * ( (lambda + 2*mu) * Dmy_2d) + ...
      (H_2d \ e_S * Hx * e_S') * (-(lambda + 2*mu) * Dmy_2d);

% Total SAT matrix
S = [S11, S12; 
     S21, S22];

A_tilde = A - S; 


% Check if  negative semi definite 

H = [H_2d, zeros(mx*my, mx*my);
    zeros(mx*my, mx*my), H_2d];

M = H * A_tilde; 


% Run test to check if A_tilde is NSD and symmetric
rel_diff = norm(M - M.', 'fro') / norm(M, 'fro');

if rel_diff < 1e-12
    disp('Matrix is symmetric');
else
    disp('Matrix is asymmetric');
    fprintf('Relative difference: %e\n', rel_diff);
end


ev = eig(full(M)); 
is_NSD = all(real(ev) <= 1e-10);

if is_NSD
    disp('Matrix is negative semidefinite');
else
    disp('Matrix is NOT negative semidefinite');
    fprintf('Largest eigenvalue: %e\n', max(real(ev)));
end



