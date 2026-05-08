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

%% --- Implementation without SAT ---

clear

% Grid parameters
mx = 80; % Grid points
my = 80;
Lx = 10; % Domain lenght
Ly = 10;
hx = Lx/(mx-1); % Step lenght 
hy = Ly/(my-1);
margin = round(mx*0.2);

x = linspace(0, Lx, mx).' ; % Define x and y vector 
y = linspace(0, Ly, my).' ;
[X, Y] = meshgrid(x, y); % Create grid in 2D

% Material constants
lambda = 22.7 * 10^9;
mu = 42.3 * 10^9;
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


%% SAT implementation

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


%% RK4 

% Time 
dt = 0.00001;          % Time step 
T_final = 0.005;       % Total time
t_steps = T_final/dt;

% Define a gaussian pulse in the centre of the domain
% Parameters
sigma = 0.9; % width
x0 = Lx/2;   % Center of x
y0 = Ly/2;   % Center of y

% Create pulse (Gaussian distribution in 2D)
pulse = exp(-( (X-x0).^2 + (Y-y0).^2 ) / (2*sigma^2));

% Place the pulse in vx and vy 
v_initial = zeros(2*mx*my, 1);
v_initial(1:mx*my) = reshape(pulse, [], 1); % For vx 
%v_initial(mx*my+1:end) = reshape(pulse, [], 1); % For vy, kolla på denna 

% Update W
w_initial = zeros(2*mx*my, 1);
u = [v_initial; w_initial];

% Definition of function f(u)
f = @(u) [u(2*mx*my+1:end); (1/rho).*A_tilde * u(1:2*mx*my)];

% --- RK4 Loop ---
for n = 1:t_steps
    k1 = f(u);
    k2 = f(u + dt/2 * k1);
    k3 = f(u + dt/2 * k2);
    k4 = f(u + dt * k3);
    
    u = u + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
    
    % Plot
    if mod(n, 1) == 0
        v_current = u(1:mx*my);
        V_plot = reshape(v_current, my, mx);
        imagesc(x, y, V_plot);
        clim([-0.5 0.5]); 
        colorbar;
        title(['Time: ', num2str(n*dt)]);
        drawnow;
    end
end