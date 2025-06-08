clear;
close all;
constants;

% 仿真时间设置
delta_t = 5e-12;
n_t     = 600;
t_start = 0;
t_end   = t_start + n_t * delta_t;
t       = t_start : delta_t : t_end;

% 初始磁矩设置
m_total = zeros(3 , length(t));
m_init  = [cos(phi_B) ; sin(phi_B) ; 0];
m_total(1:3 , 1) = m_init;

% 四阶Runge-Kutta方法求解
for i = 1:n_t
    
    t_curr = t(i);
    
    k1 = Dmdt(m_total(:,i) , t_curr);
    k2 = Dmdt(m_total(:,i) + delta_t * k1 / 2 , t_curr + delta_t / 2);
    k3 = Dmdt(m_total(:,i) + delta_t * k2 / 2 , t_curr + delta_t / 2);
    k4 = Dmdt(m_total(:,i) + delta_t * k3 , t_curr + delta_t);
    
    delta_m_total = delta_t * (k1 + 2*k2 + 2*k3 + k4) / 6;
    
    m_total(: , i+1) = m_total(:,i) + delta_m_total;
    
    % 每次循环对磁矩方向归一化，防止数值耗散
    m_total(: , i+1) = m_total(: , i+1) / norm(m_total(: , i+1));
    
end

% 坐标变换
m_total_trans = [cos(phi_B) sin(phi_B) 0 ; -sin(phi_B) cos(phi_B) 0 ; 0 0 1] * m_total;

omega = 2 * pi * f;
B_OF  = B_Oe + B_FL;

% 求解线性方程组
A = [ omega*1i , B_ext*gamma + alpha*omega*1i + M_s*gamma*mu_0; -B_ext*gamma - alpha*omega*1i , omega*1i];
b = [-B_DL*gamma*cos(phi_B) ; B_OF*gamma*cos(phi_B)];
m = A^-1 * b;
my = m(1);
mz = m(2);

% my = -(gamma*cos(phi_B)*(B_DL*omega*1i + B_OF*B_ext*gamma + B_OF*alpha*omega*1i + B_OF*M_s*gamma*mu_0))/...
%     (B_ext^2*gamma^2 + B_ext*alpha*gamma*omega*2i + M_s*mu_0*B_ext*gamma^2 - alpha^2*omega^2 + M_s*mu_0*alpha*gamma*omega*1i - omega^2);
% mz = -(gamma*cos(phi_B)*(B_DL*B_ext*gamma - B_OF*omega*1i + B_DL*alpha*omega*1i))/...
%     (B_ext^2*gamma^2 + B_ext*alpha*gamma*omega*2i + M_s*mu_0*B_ext*gamma^2 - alpha^2*omega^2 + M_s*mu_0*alpha*gamma*omega*1i - omega^2);

figure;
subplot(311)
plot(t,m_total_trans(1,:),t,sqrt(1 - real(my * exp(1i * omega * t)).^2 - real(mz * exp(1i * omega * t)).^2))

subplot(312)
plot(t,m_total_trans(2,:),t,real(my * exp(1i * omega * t)))

subplot(313)
plot(t,m_total_trans(3,:),t,real(mz * exp(1i * omega * t)))