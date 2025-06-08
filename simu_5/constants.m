global B_x...
    B_y   ...
    B_z   ...
    N_x   ...
    N_y   ...
    N_z   ...
    mu_0  ...
    M_s   ;

% 设置磁场强度和方向
B_ext   = 0.137;
phi_B = 0.2;

B_x     = B_ext * cos(phi_B);
B_y     = B_ext * sin(phi_B);
B_z     = 0;
N_x     = 0;
N_y     = 0;
N_z     = 1;
mu_0    = 4 * pi * 1e-7;
M_s     = 6.4e+5;


global gamma...
    sigma   ...
    y       ...
    alpha   ;

u_B      = 9.274e-24;
hbar     = 6.62607015e-34 / (2 * pi);
gamma    = 2 * u_B / hbar;
sigma    = [0 ; -1 ; 0];
y        = [0 ; 1 ; 0];
alpha    = 0.02;


global f ...
    B_Oe ...
    B_FL ...
    B_DL ;

f     = 5e+9;
B_Oe  = 4 * 1e-4;
B_FL  = 1 * 1e-4;
B_DL  = 10 * 1e-4;