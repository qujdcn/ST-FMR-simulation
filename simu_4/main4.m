clear;
close all;

% 模拟参数
mu_0  = 4 * pi * 1e-7;
M_s   = 6.4e+5;

u_B   = 9.274e-24;
hbar  = 6.62607015e-34 / (2 * pi);
gamma = 2 * u_B / hbar;
alpha = 0.02;

f     = 10e+9;
omega = 2 * pi * f;

B_Oe  = 4 * 1e-4;
B_FL  = 1 * 1e-4;
B_DL  = 10 * 1e-4;
B_OF  = B_Oe + B_FL;

phi_B = 0.1;
B_ext = 0.00:0.001:0.3;
n     = length(B_ext);
mx1   = zeros(1 , n);

% 扫描外加磁场Bext
for i = 1:n
    
    B_x = B_ext(i) * cos(phi_B);
    B_y = B_ext(i) * sin(phi_B);
    
    % 计算线性方程组，求解磁矩的进动状态
    A = [1i * omega , 0 , -gamma * B_y - mu_0 * M_s * gamma * sin(phi_B) - 1i * omega * alpha * sin(phi_B);
        0 , 1i * omega , gamma * B_x + mu_0 * M_s * gamma * cos(phi_B) + 1i * omega * alpha * cos(phi_B);
        gamma * B_y + 1i * omega * alpha * cos(phi_B) , -gamma * B_x - 1i * omega * alpha * sin(phi_B) , 1i * omega];
    b = [gamma * B_DL * cos(phi_B) * sin(phi_B) ; -gamma * B_DL * cos(phi_B)^2 ; gamma * B_OF * cos(phi_B)];
    
    m = A^-1 * b;
    mx1(i) = m(1);
    
end

% 计算共振场（求根公式）
B_r = (-gamma^2 * mu_0 * M_s + sqrt((gamma^2 * mu_0 * M_s)^2 + 4 * gamma^2 * omega^2)) / (2 * gamma^2);
% 利用公式求解real(mx1)
real_mx1 = sin(phi_B) * cos(phi_B) * (omega^2 * alpha * B_DL + gamma^2 * (B_r + mu_0 * M_s) * (B_ext - B_r) * B_OF) ./...
    ((2 * B_r + mu_0 * M_s) * (omega^2 * alpha^2 + gamma^2 * (B_ext - B_r).^2));

% 对比两种计算方法的计算结果
figure;
hold on
plot(B_ext , real(mx1),'k')
plot(B_ext , real_mx1,'r')