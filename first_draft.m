%% ---Here we have all our functions---

function [Dp, Dm, H] = SBP_5(m, h) 
    H=diag(ones(m,1),0);
    H(1:4,1:4)=[0.251e3 / 0.720e3 0 0 0; 0 0.299e3 / 0.240e3 0 0; 0 0 0.211e3 / 0.240e3 0; 0 0 0 0.739e3 / 0.720e3;];
    H(m-3:m,m-3:m)=fliplr(flipud(H(1:4,1:4)));
    H=H*h;
    HI=inv(H);
    
    Qp=(1/20*diag(ones(m-2,1),-2)-1/2*diag(ones(m-1,1),-1)-1/3*diag(ones(m,1),0)+1*diag(ones(m-1,1),+1)-1/4*diag(ones(m-2,1),+2)+1/30*diag(ones(m-3,1),+3));

    Q_U = [-0.1e1 / 0.120e3 0.941e3 / 0.1440e4 -0.47e2 / 0.360e3 -0.7e1 / 0.480e3; -0.869e3 / 0.1440e4 -0.11e2 / 0.120e3 0.25e2 / 0.32e2 -0.43e2 / 0.360e3; 0.29e2 / 0.360e3 -0.17e2 / 0.32e2 -0.29e2 / 0.120e3 0.1309e4 / 0.1440e4; 0.1e1 / 0.32e2 -0.11e2 / 0.360e3 -0.661e3 / 0.1440e4 -0.13e2 / 0.40e2;];

    Qp(1:4,1:4)=Q_U;
    Qp(m-3:m,m-3:m)=flipud( fliplr(Q_U(1:4,1:4) ) )'; %%% This is different from standard SBP

    Qm=-Qp';

    e_1=zeros(m,1);e_1(1)=1;
    e_m=zeros(m,1);e_m(m)=1;

    Dp=HI*(Qp-1/2*e_1*e_1'+1/2*e_m*e_m') ;

    Dm=HI*(Qm-1/2*e_1*e_1'+1/2*e_m*e_m') ;
end

function [Dp, Dm, H] = SBP_7(m, h) 
    H=diag(ones(m,1),0);
    H(1:6,1:6)=[0.19087e5 / 0.60480e5 0 0 0 0 0; 0 0.84199e5 / 0.60480e5 0 0 0 0; 0 0 0.18869e5 / 0.30240e5 0 0 0; 0 0 0 0.37621e5 / 0.30240e5 0 0; 0 0 0 0 0.55031e5 / 0.60480e5 0; 0 0 0 0 0 0.61343e5 / 0.60480e5;];
    H(m-5:m,m-5:m)=fliplr(flipud(H(1:6,1:6)));
    H=H*h;
    HI=inv(H);

    Qp=(-1/105*diag(ones(m-3,1),-3)+1/10*diag(ones(m-2,1),-2)-3/5*diag(ones(m-1,1),-1)-1/4*diag(ones(m,1),0)+1*diag(ones(m-1,1),+1)-3/10*diag(ones(m-2,1),+2)+1/15*diag(ones(m-3,1),+3)-1/140*diag(ones(m-4,1),+4));
    Q_U =[-0.265e3 / 0.300272e6 0.1587945773e10 / 0.2432203200e10 -0.1926361e7 / 0.25737600e8 -0.84398989e8 / 0.810734400e9 0.48781961e8 / 0.4864406400e10 0.3429119e7 / 0.202683600e9; -0.1570125773e10 / 0.2432203200e10 -0.26517e5 / 0.1501360e7 0.240029831e9 / 0.486440640e9 0.202934303e9 / 0.972881280e9 0.118207e6 / 0.13512240e8 -0.231357719e9 / 0.4864406400e10; 0.1626361e7 / 0.25737600e8 -0.206937767e9 / 0.486440640e9 -0.61067e5 / 0.750680e6 0.49602727e8 / 0.81073440e8 -0.43783933e8 / 0.194576256e9 0.51815011e8 / 0.810734400e9; 0.91418989e8 / 0.810734400e9 -0.53314099e8 / 0.194576256e9 -0.33094279e8 / 0.81073440e8 -0.18269e5 / 0.107240e6 0.440626231e9 / 0.486440640e9 -0.365711063e9 / 0.1621468800e10; -0.62551961e8 / 0.4864406400e10 0.799e3 / 0.35280e5 0.82588241e8 / 0.972881280e9 -0.279245719e9 / 0.486440640e9 -0.346583e6 / 0.1501360e7 0.2312302333e10 / 0.2432203200e10; -0.3375119e7 / 0.202683600e9 0.202087559e9 / 0.4864406400e10 -0.11297731e8 / 0.810734400e9 0.61008503e8 / 0.1621468800e10 -0.1360092253e10 / 0.2432203200e10 -0.10677e5 / 0.42896e5;];

    Qp(1:6,1:6)=Q_U;
    Qp(m-5:m,m-5:m)=flipud( fliplr(Q_U(1:6,1:6) ) )'; %%% This is different from standard SBP

    Qm=-Qp';

    e_1=zeros(m,1);e_1(1)=1;
    e_m=zeros(m,1);e_m(m)=1;

    Dp=HI*(Qp-1/2*e_1*e_1'+1/2*e_m*e_m') ;

    Dm=HI*(Qm-1/2*e_1*e_1'+1/2*e_m*e_m') ;
end

function [Dp, Dm, H] = SBP_9(m, h) 
    H=diag(ones(m,1),0);
    H(1:8,1:8)=[0.1070017e7 / 0.3628800e7 0 0 0 0 0 0 0; 0 0.5537111e7 / 0.3628800e7 0 0 0 0 0 0; 0 0 0.103613e6 / 0.403200e6 0 0 0 0 0; 0 0 0 0.261115e6 / 0.145152e6 0 0 0 0; 0 0 0 0 0.298951e6 / 0.725760e6 0 0 0; 0 0 0 0 0 0.515677e6 / 0.403200e6 0 0; 0 0 0 0 0 0 0.3349879e7 / 0.3628800e7 0; 0 0 0 0 0 0 0 0.3662753e7 / 0.3628800e7;];

    H(m-7:m,m-7:m)=fliplr(flipud(H(1:8,1:8)));
    H=H*h;
    HI=inv(H);
    
    Qp=(1/504*diag(ones(m-4,1),-4)-1/42*diag(ones(m-3,1),-3)+1/7*diag(ones(m-2,1),-2)-2/3*diag(ones(m-1,1),-1)-1/5*diag(ones(m,1),0)+1*diag(ones(m-1,1),+1)-1/3*diag(ones(m-2,1),+2)+2/21*diag(ones(m-3,1),+3)-1/56*diag(ones(m-4,1),+4)+1/630*diag(ones(m-5,1),+5));
    Q_U =[-0.5561e4 / 0.47263920e8 0.4186300102421e13 / 0.6193464076800e13 -0.377895002003e12 / 0.5806372572000e13 -0.16485548951749e14 / 0.111482353382400e15 -0.113245973003e12 / 0.3716078446080e13 0.355360297339e12 / 0.4645098057600e13 0.321012170669e12 / 0.55741176691200e14 -0.388397049437e12 / 0.26543417472000e14; -0.4178798062421e13 / 0.6193464076800e13 -0.493793e6 / 0.141791760e9 0.725405227507e12 / 0.2413037952000e13 0.3904159533697e13 / 0.9290196115200e13 0.2483046570341e13 / 0.13935294172800e14 -0.4336328670953e13 / 0.18580392230400e14 -0.1258688487061e13 / 0.37160784460800e14 0.12931584852209e14 / 0.278705883456000e15; 0.363359390003e12 / 0.5806372572000e13 -0.7539548734577e13 / 0.26543417472000e14 -0.69332623e8 / 0.2977626960e10 0.9994352248429e13 / 0.18580392230400e14 -0.8195655811631e13 / 0.18580392230400e14 0.7361486640463e13 / 0.61934640768000e14 0.5539855071347e13 / 0.92901961152000e14 -0.12898722943e11 / 0.422281641600e12; 0.16773595838149e14 / 0.111482353382400e15 -0.372477950627e12 / 0.844563283200e12 -0.8659050093229e13 / 0.18580392230400e14 -0.207799621e9 / 0.2977626960e10 0.1734921317461e13 / 0.2477385630720e13 0.2530020015841e13 / 0.18580392230400e14 0.441856623253e12 / 0.13935294172800e14 -0.115132773073e12 / 0.2654341747200e13; 0.108449122763e12 / 0.3716078446080e13 -0.2283566671541e13 / 0.13935294172800e14 0.6976424333231e13 / 0.18580392230400e14 -0.440819477447e12 / 0.825795210240e12 -0.55386253e8 / 0.425375280e9 0.2479572560009e13 / 0.3716078446080e13 -0.40258468963e11 / 0.120651897600e12 0.11808221047099e14 / 0.111482353382400e15; -0.32231128289e11 / 0.422281641600e12 0.4244793299753e13 / 0.18580392230400e14 -0.5173673584463e13 / 0.61934640768000e14 -0.4848139955041e13 / 0.18580392230400e14 -0.1506045711689e13 / 0.3716078446080e13 -0.526653889e9 / 0.2977626960e10 0.36411368691307e14 / 0.37160784460800e14 -0.825434105779e12 / 0.2903186286000e13; -0.316459841069e12 / 0.55741176691200e14 0.1277069729941e13 / 0.37160784460800e14 -0.6499182375347e13 / 0.92901961152000e14 0.355606625147e12 / 0.13935294172800e14 0.1519272420551e13 / 0.9290196115200e13 -0.2240079855137e13 / 0.3378253132800e13 -0.584765899e9 / 0.2977626960e10 0.2301241355533e13 / 0.2382101568000e13; 0.387779289437e12 / 0.26543417472000e14 -0.12908508708209e14 / 0.278705883456000e15 0.147710908133e12 / 0.4645098057600e13 0.534025841911e12 / 0.18580392230400e14 -0.4119981443899e13 / 0.111482353382400e15 0.279819152779e12 / 0.2903186286000e13 -0.1510324515533e13 / 0.2382101568000e13 -0.85017967e8 / 0.425375280e9;];

    Qp(1:8,1:8)=Q_U;
    Qp(m-7:m,m-7:m)=flipud( fliplr(Q_U(1:8,1:8) ) )'; %%% This is different from standard SBP

    Qm=-Qp';

    e_1=zeros(m,1);e_1(1)=1;
    e_m=zeros(m,1);e_m(m)=1;

    Dp=HI*(Qp-1/2*e_1*e_1'+1/2*e_m*e_m') ;

    Dm=HI*(Qm-1/2*e_1*e_1'+1/2*e_m*e_m') ;
end

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
