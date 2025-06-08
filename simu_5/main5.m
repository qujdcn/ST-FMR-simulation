clear;

syms my mz phi_B omega B_ext gamma alpha B_OF B_DL mu_0 M_s;
y = [sin(phi_B) ; cos(phi_B) ; 0];
sigma = -y;
m = [1 ; my ; mz];


dmdt = [0 ; 1i * omega * my ; 1i * omega * mz];

% cross(m , [0 ; 0 ; 1] .* m)
right1 = gamma * mu_0 * M_s * [0 ; -mz ; 0] - gamma * cross(m,[B_ext ; 0 ; 0]);

right2 = alpha * cross(m , dmdt);

% cross(sigma , m)
right34 = gamma * B_OF * [0 ; 0 ; cos(phi_B)];

% cross(m , cross(sigma , m))
right5 = gamma * B_DL * [0 ; -cos(phi_B) ; 0];

ALL = dmdt - (right1 + right2 + right34 + right5);

A = [ omega*1i , B_ext*gamma + alpha*omega*1i + M_s*gamma*mu_0; -B_ext*gamma - alpha*omega*1i , omega*1i];
b = [-B_DL*gamma*cos(phi_B) ; B_OF*gamma*cos(phi_B)];

simplify(A ^ -1 * b)