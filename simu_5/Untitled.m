clear;

syms my mz phi_B omega B_ext gamma alpha B_OF B_DL mu_0 M_s;
y = [sin(phi_B) ; cos(phi_B) ; 0];
sigma = -y;
m = [1 ; my ; mz];

% cross(sigma , m)
% cross(m , cross(sigma , m))

dmdt = [0 ; 1i * omega * my ; 1i * omega * mz];

left2 = alpha * cross(m , dmdt);

left34 = gamma * B_OF * [0 ; 0 ; cos(phi_B)];

left5 = gamma * B_DL * [0 ; -cos(phi_B) ; 0];

left1 = gamma * mu_0 * M_s * [0 ; -mz ; 0] - gamma * cross(m,[B_ext ; 0 ; 0]);

ALL = left1 + left2 + left34 + left5 - dmdt;



% my = (B_OF*gamma*cos(phi))/(B*gamma + alpha*omega*1i + gamma*muMs) + (B_DL*gamma*cos(phi)*1i)/omega;
% mz = -(B_OF*gamma*cos(phi)*1i)/omega + (B_DL*gamma*cos(phi))/(B*gamma + alpha*omega*1i);
A = [ -omega*1i , -B_ext*gamma - alpha*omega*1i - M_s*gamma*mu_0; B_ext*gamma + alpha*omega*1i , -omega*1i];
b = [B_DL*gamma*cos(phi_B) ; -B_OF*gamma*cos(phi_B)];
% 
A ^ -1 * b