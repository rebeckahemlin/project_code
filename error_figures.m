%% --- Implementation without SAT ---

clear

% Grid parameters
mx = 100; % Grid points
my = 100;
Lx = 10; % Domain lenght
Ly = 10;
hx = Lx/(mx-1); % Step lenght 
hy = Ly/(my-1);
margin = round(mx*0.2);

x = linspace(0, Lx, mx).' ; % Define x and y vector 
y = linspace(0, Ly, my).' ;
[X, Y] = meshgrid(x, y); % Create grid i 2D

% Material constants
lambda = 42.3* 10^9;
mu = 22.7 * 10^9;
roh = 2700;

% Create 1D operators for x and y 
[Dpx, Dmx, Hx] = SBP_9(mx, hx);
[Dpy, Dmy, Hy] = SBP_9(my, hy);

% Create identity matrices 
Ix = eye(mx); 
Iy = eye(my);

Dpx_2d = kron(sparse(Dpx), sparse(Iy)); % Derivative in x direction with Dp
Dpy_2d = kron(sparse(Ix), sparse(Dpy)); % Derivative in y direction wuth Dp

Dmx_2d = kron(sparse(Dmx), sparse(Iy)); % Derivative in x direction with Dm
Dmy_2d = kron(sparse(Ix), sparse(Dmy)); % Derivative in x direction with Dm

% Elements in derivative matrix A
A11 = Dpx_2d * (lambda + 2*mu) * Dmx_2d + Dpy_2d * mu * Dmy_2d;
A12 = Dpy_2d * mu * Dmx_2d + Dpx_2d * lambda * Dmy_2d; 
A21 = Dpx_2d * mu * Dmy_2d + Dpy_2d * lambda * Dmx_2d;
A22 = Dpx_2d * mu * Dmx_2d + Dpy_2d * (lambda + 2*mu) * Dmy_2d;

A = [A11, A12;
    A21, A22];

% Analytical solution 
syms x y lambda mu
syms vx(x,y) vy(x,y)

% Symbolic differentiation 
PDE_x_symbolic = diff((lambda + 2*mu)*diff(vx, x), x) + diff(mu*diff(vy, x), y) + diff(mu*diff(vx, y), y) + diff(lambda*diff(vy, y), x);
PDE_y_symbolic = diff((lambda + 2*mu)*diff(vy, y), y) + diff(mu*diff(vx, y), x) + diff(mu*diff(vy, x), x) + diff(lambda*diff(vx, x), y);

% test functions
vx_test = sin(x) * sin(y); 
vy_test = cos(x) * cos(y); 
V_test = [vx_test; vy_test];

% Analytical solutions to our test functions
PDE_x_analytical = subs(PDE_x_symbolic, {vx(x,y), vy(x,y)}, {vx_test, vy_test});
PDE_y_analytical = subs(PDE_y_symbolic, {vx(x,y), vy(x,y)}, {vx_test, vy_test});

PDE_x_analytical_function = matlabFunction(PDE_x_analytical, 'Vars', [x, y, lambda, mu]);
PDE_y_analytical_function = matlabFunction(PDE_y_analytical, 'Vars', [x, y, lambda, mu]);

% Now we put in our dicrete X and y and our value for lambda and mu
analytical_result_x = PDE_x_analytical_function(X, Y, 42.3 * 10^9,  22.7 * 10^9); % x, y, lambda, mu
analytical_result_y = PDE_y_analytical_function(X, Y, 42.3 * 10^9, 22.7 * 10^9);

% Now solve for the numerical solution of test function with A
vx_discrete = sin(X) .* sin(Y);
vy_discrete = cos(X) .* cos(Y);

%flatten matrices into vectors
vx_vec = vx_discrete(:);
vy_vec = vy_discrete(:);

V_num = [vx_vec; vy_vec];

% calculate RHS
result_num = A * V_num;

% shape back into matrices
numerical_result_x = reshape(result_num(1:mx*my), my, mx);
numerical_result_y = reshape(result_num(mx*my+1:end), my, mx);

% calculate difference between analytical and numerical
error_x = abs(analytical_result_x - numerical_result_x);
error_y = abs(analytical_result_y - numerical_result_y);

fprintf('Maximum error in x: %e\n', max(error_x(:)));
fprintf('Maximum error in y: %e\n', max(error_y(:)));


% Visulazation 
figure('Name', '2D SBP Test', 'Position', [100, 100, 1200, 400]);

subplot(1,2,1);
surf(X, Y, error_x);
set(gca, 'ZScale', 'log'); % Logaritmic scale 
title('Absolute error in df/dx (Logaritmic scale)');
xlabel('x'); ylabel('y'); shading interp; colorbar;

subplot(1,2,2);
surf(X, Y, error_y);
set(gca, 'ZScale', 'log'); % Logaritmic scale 
title('Absolute error in df/dy (Logaritmic scale)');
xlabel('x'); ylabel('y'); shading interp; colorbar; 

