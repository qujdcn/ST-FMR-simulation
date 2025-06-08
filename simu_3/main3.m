clear;
close all;
constants;

% ����ʱ������
delta_t = 5e-12;
n_t     = 600;
t_start = 0;
t_end   = t_start + n_t * delta_t;
t       = t_start : delta_t : t_end;

% ��ʼ�ž�����
m_total = zeros(3 , length(t));
m_init  = [cos(phi_B) ; sin(phi_B) ; 0];
m_total(1:3 , 1) = m_init;

% �Ľ�Runge-Kutta�������
for i = 1:n_t
    
    t_curr = t(i);
    
    k1 = Dmdt(m_total(:,i) , t_curr);
    k2 = Dmdt(m_total(:,i) + delta_t * k1 / 2 , t_curr + delta_t / 2);
    k3 = Dmdt(m_total(:,i) + delta_t * k2 / 2 , t_curr + delta_t / 2);
    k4 = Dmdt(m_total(:,i) + delta_t * k3 , t_curr + delta_t);
    
    delta_m_total = delta_t * (k1 + 2*k2 + 2*k3 + k4) / 6;
    
    m_total(: , i+1) = m_total(:,i) + delta_m_total;
    
    % ÿ��ѭ���Դžط����һ������ֹ��ֵ��ɢ
    m_total(: , i+1) = m_total(: , i+1) / norm(m_total(: , i+1));
    
end

% �žؽ�����
mx0 = cos(phi_B);
my0 = sin(phi_B);
mz0 = 0;

% �������Է����飬���žصĽ���״̬
omega = 2 * pi * f;
B_OF = B_Oe + B_FL;

A = [1i * omega , 0 , -gamma * B_y - mu_0 * M_s * gamma * sin(phi_B) - 1i * omega * alpha * sin(phi_B);
    0 , 1i * omega , gamma * B_x + mu_0 * M_s * gamma * cos(phi_B) + 1i * omega * alpha * cos(phi_B);
    gamma * B_y + 1i * omega * alpha * cos(phi_B) , -gamma * B_x - 1i * omega * alpha * sin(phi_B) , 1i * omega];
b = [gamma * B_DL * cos(phi_B) * sin(phi_B) ; -gamma * B_DL * cos(phi_B)^2 ; gamma * B_OF * cos(phi_B)];

m = A^-1 * b;
mx1 = m(1);
my1 = m(2);
mz1 = -(gamma*cos(phi_B)*(B_DL*B_ext*gamma - B_OF*omega*1i + B_DL*alpha*omega*1i))/...
    (B_ext^2*gamma^2 + B_ext*alpha*gamma*omega*2i + M_s*mu_0*B_ext*gamma^2 - alpha^2*omega^2 + M_s*mu_0*alpha*gamma*omega*1i - omega^2);

% ����mx��my��mz��ʱ��ı仯ͼ�񣬶Ա����ּ��㷽���ļ�����
figure;
subplot(311)
plot(t,m_total(1,:) , t,real(mx0 + mx1 * exp(1i * omega * t)))

subplot(312)
plot(t,m_total(2,:) , t,real(my0 + my1 * exp(1i * omega * t)))

subplot(313)
plot(t,m_total(3,:) , t,real(mz0 + mz1 * exp(1i * omega * t)))