%% --- Testing of SPB operators ---

mx = 100; % Grid points
my = 100;
Lx = 10; % Domain lenght
Ly = 10;
hx = Lx/(mx-1); % Step lenght 
hy = Ly/(my-1);

margin = round(mx*0.2);

x = linspace(0, Lx, mx).' ; % Define x and y vector 
y = linspace(0, Ly, my).' ;
[X, Y] = meshgrid(x, y); % Create grid in 2D

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

% Functions to test
F = sin(X) .* e^Y + 5;
dFdx_exact = cos(X) .* e^Y;
dFdy_exact = sin(X) .* e^Y;

% Numerical solutions of our test function
f_vec = F(:); % Convert my x mx matrix to  my*mx x 1 vector

%% ---Numerical differentiation solution with Dp---
dfdy_vec = Dpy_2d * f_vec; 
dfdx_vec = Dpx_2d * f_vec;

DFDY_num = reshape(dfdy_vec, my, mx); % Convert back to my x mx matrix 
DFDX_num = reshape(dfdx_vec, my, mx);

% Error calculations 
error_x = abs(DFDX_num - dFdx_exact);
error_y = abs(DFDY_num - dFdy_exact);


% Visulazation 
figure('Name', '2D SBP Test', 'Position', [100, 100, 1200, 400]);

subplot(1,3,1);
surf(X, Y, F);
title('f(x,y) = sin(x)cos(y)');
xlabel('x'); ylabel('y'); shading interp; colorbar;

subplot(1,3,2);
surf(X, Y, error_x);
set(gca, 'ZScale', 'log'); % Logaritmic scale 
title('Absolute error (Logaritmic scale)');
xlabel('x'); ylabel('y'); shading interp; colorbar;

subplot(1,3,3);
surf(X, Y, error_y);
set(gca, 'ZScale', 'log'); % Logaritmic scale 
title('Absolute error (Logaritmic scale)');
xlabel('x'); ylabel('y'); shading interp; colorbar; 


fprintf('Maximum error in x derivative with Dp: %e\n', max(error_x(:)))
fprintf('Maximum error in y derivative with Dp: %e\n', max(error_y(:)))

inner_error_x = error_x(margin:end-margin, margin:end-margin);
inner_error_y = error_y(margin:end-margin, margin:end-margin);

rms_inner_x = sqrt(mean(inner_error_x(:).^2));
rms_inner_y = sqrt(mean(inner_error_y(:).^2));

fprintf('RMS Error in domain center, in x derivative with Dp: %e\n', rms_inner_x)
fprintf('RMS Error in domain center, in y derivative with Dp: %e\n', rms_inner_y)


%% ---Numerical differentiation solution with Dm---
dfdy_vec = Dmy_2d * f_vec; 
dfdx_vec = Dmx_2d * f_vec;

DFDY_num = reshape(dfdy_vec, my, mx); % Convert back to a my x mx matrix 
DFDX_num = reshape(dfdx_vec, my, mx);

% Error calculations 
error_x = abs(DFDX_num - dFdx_exact);
error_y = abs(DFDY_num - dFdy_exact);


% Visulazation 
figure('Name', '2D SBP Test', 'Position', [100, 100, 1200, 400]);

subplot(1,3,1);
surf(X, Y, F);
title('f(x,y) = sin(x)cos(y)');
xlabel('x'); ylabel('y'); shading interp; colorbar;

subplot(1,3,2);
surf(X, Y, error_x);
set(gca, 'ZScale', 'log'); % Logaritmic scale 
title('Absolute error (Logaritmic scale)');
xlabel('x'); ylabel('y'); shading interp; colorbar;

subplot(1,3,3);
surf(X, Y, error_y);
set(gca, 'ZScale', 'log'); % Logaritmic scale
title('Absolute error in (Logaritmic scale)');
xlabel('x'); ylabel('y'); shading interp; colorbar;


fprintf('Maximum error in x derivative with Dm: %e\n', max(error_x(:)))
fprintf('Maximum error in y derivative with Dm: %e\n', max(error_y(:)))

inner_error_x = error_x(margin:end-margin, margin:end-margin);
inner_error_y = error_y(margin:end-margin, margin:end-margin);

rms_inner_x = sqrt(mean(inner_error_x(:).^2));
rms_inner_y = sqrt(mean(inner_error_y(:).^2));

fprintf('RMS Error in domain center, in x derivative with Dm: %e\n', rms_inner_x)
fprintf('RMS Error in domain center, in y derivative with Dm: %e\n', rms_inner_y)
%% Plotting error

% list with error

m_list = [12, 25, 50, 100];
h_list = 10 ./ (m_list - 1);


max_err_list_Dp_x_7 = [3.215983e-01, 2.487731e-02, 1.597694e-03, 2.219221e-04];
max_err_list_Dp_y_7 = [1.805455e-01, 3.761110e-02, 4.540804e-03, 5.561048e-04];
max_err_list_Dm_x_7 = [3.249224e-01, 2.511192e-02, 1.612337e-03, 2.200576e-04];
max_err_list_Dm_y_7 = [1.790851e-01, 3.727808e-02, 4.501379e-03, 5.609116e-04];

max_err_list_Dp_x_9 = [3.082654e-01, 2.049917e-02, 1.470424e-03, 8.635498e-05];
max_err_list_Dp_y_9 = [5.722413e-01, 2.140558e-02, 7.653787e-04, 2.683994e-05];
max_err_list_Dm_x_9 = [3.082654e-01, 2.048124e-02, 1.469383e-03, 8.641642e-05];
max_err_list_Dm_y_9 = [5.722413e-01, 2.141999e-02, 7.658084e-04, 2.681733e-05];

centrum_err_list_Dp_x_9 = [1.069679e-01, 7.949623e-04, 2.404549e-10, 4.327732e-13];
centrum_err_list_Dp_y_9 = [1.595148e-01, 1.526423e-03, 2.405924e-10, 4.323654e-13];
centrum_err_list_Dm_x_9 = [1.062950e-01, 7.916315e-04, 2.405652e-10, 4.323880e-13];
centrum_err_list_Dm_y_9 = [1.827908e-01, 1.526608e-03, 2.404807e-10, 4.327453e-13];



%% --- Plotting for maximum error ---
ref_h3 = max_err_list_Dp_x_9(1) * (h_list / h_list(1)).^4;
ref_h4 = max_err_list_Dp_x_9(1) * (h_list / h_list(1)).^5;

figure;             
loglog(h_list, max_err_list_Dp_x_9, '-or', 'LineWidth', 1.5); hold on;
loglog(h_list, max_err_list_Dp_y_9, '-sq', 'LineWidth', 1.5);
loglog(h_list, max_err_list_Dm_x_9, '-^b', 'LineWidth', 1.5);
loglog(h_list, max_err_list_Dm_y_9, '-*m', 'LineWidth', 1.5);

loglog(h_list, ref_h3, '--k', 'LineWidth', 2)
loglog(h_list, ref_h4, ':k', 'LineWidth', 1.5, 'DisplayName', 'O(h^4) Slope');

grid on;
xlabel('Step h');
ylabel('Maximal error');
title('Convergence analysis');
legend('D+ x', 'D+ y', 'D- x', 'D- y', 'Location', 'best');

%% --- Plotting for centrum error ---
ref_h12= centrum_err_list_Dp_x_9(1) * (h_list / h_list(1)).^9;
ref_h11= centrum_err_list_Dp_x_9(1) * (h_list / h_list(1)).^8;

figure;             
loglog(h_list, centrum_err_list_Dp_x_9, '-or', 'LineWidth', 1.5); hold on;
loglog(h_list, centrum_err_list_Dp_y_9, '-sq', 'LineWidth', 1.5);
loglog(h_list, centrum_err_list_Dm_x_9, '-^b', 'LineWidth', 1.5);
loglog(h_list, centrum_err_list_Dm_y_9, '-*m', 'LineWidth', 1.5);

loglog(h_list, ref_h12, '--k', 'LineWidth', 2)
loglog(h_list, ref_h11, ':k', 'LineWidth', 1.5, 'DisplayName', 'O(h^4) Slope');

grid on;
xlabel('Steglängd h');
ylabel('Maximalt fel');
title('Konvergensanalys för olika operatorer');
legend('D+ x', 'D+ y', 'D- x', 'D- y', 'Location', 'best');
