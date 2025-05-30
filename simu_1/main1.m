clear;
close all;
constants;

% ����ʱ������
delta_t = 5e-12;
n_t     = 2000;
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

% ����mx��my��mz��ʱ��ı仯ͼ��
figure;
subplot(311)
plot(t,m_total(1,:))

subplot(312)
plot(t,m_total(2,:))

subplot(313)
plot(t,m_total(3,:))