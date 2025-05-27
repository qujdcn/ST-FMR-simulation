clear;
close all;
constants;

delta_t = 5e-12;
n_t     = 2000;
t_start = 0;
t_end   = t_start + n_t * delta_t;
t       = t_start : delta_t : t_end;

m_total = zeros(3 , length(t));
m_init  = [cos(theta_B) ; sin(theta_B) ; 0];
m_total(1:3 , 1) = m_init;

for i = 1:n_t
    
    t_curr = t(i);
    
    k1 = m_prime(m_total(:,i) , t_curr);
    k2 = m_prime(m_total(:,i) + delta_t * k1 / 2 , t_curr + delta_t / 2);
    k3 = m_prime(m_total(:,i) + delta_t * k2 / 2 , t_curr + delta_t / 2);
    k4 = m_prime(m_total(:,i) + delta_t * k3 , t_curr + delta_t);
    
    delta_m_total = delta_t * (k1 + 2*k2 + 2*k3 + k4) / 6;
    
    m_total(: , i+1) = m_total(:,i) + delta_m_total;
    
    m_total(: , i+1) = m_total(: , i+1) / norm(m_total(: , i+1));
    
end


figure;
subplot(311)
plot(t,m_total(1,:),t,cos(theta_B) + real(exp(1i * 2 * pi * f *t) * m(1)))

subplot(312)
plot(t,m_total(2,:),t,sin(theta_B) + real(exp(1i * 2 * pi * f *t) * m(2)))

subplot(313)
plot(t,m_total(3,:),t,real(exp(1i * 2 * pi * f *t) * m(3)))
