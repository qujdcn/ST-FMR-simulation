% clear;
% load('matlab.mat');
% 
% dmdt = zeros(size(m_total));
% 
% dmdt(: , 2:600) = (m_total(: , 3:601) - m_total(: , 1:599)) / (2 * delta_t);
% 
% left1  = zeros(size(m_total));
% left2  = zeros(size(m_total));
% left34 = zeros(size(m_total));
% left5  = zeros(size(m_total));
% 
% for i = 1:601
%    
%     left1(: , i)  = -gamma * cross(m_total(: , i) , ([B_ext * cos(phi_B) ; B_ext * sin(phi_B) ; 0] + [0 ; 0 ; -mu_0 * M_s] .* m_total(: , i)));
%     left2(: , i)  = alpha * cross(m_total(: , i) , dmdt(: , i));
%     left34(: , i) = gamma * B_OF * cos(2 * pi * f * t(i)) * cross([0 ; -1 ; 0] , m_total(: , i));
%     left5(: , i)  = gamma * B_DL * cos(2 * pi * f * t(i)) * cross(m_total(: , i) , cross([0 ; -1 ; 0] , m_total(: , i)));
%     
% end

clear;
load('matlab.mat');
load('matlab2.mat');
omega = 2 * pi * f;

m_total = [cos(phi_B) sin(phi_B) 0 ; -sin(phi_B) cos(phi_B) 0 ; 0 0 1] * m_total;
figure;
subplot(311)
plot(t,m_total(1,:))

subplot(312)
plot(t,m_total(2,:),t,real(my * exp(1i * 2 * pi * f *t)))

subplot(313)
plot(t,m_total(3,:),t,real(mz * exp(1i * 2 * pi * f *t)))

dmdt = zeros(size(m_total));

dmdt(: , 2:600) = (m_total(: , 3:601) - m_total(: , 1:599)) / (2 * delta_t);

left1  = zeros(size(m_total));
left2  = zeros(size(m_total));
left34 = zeros(size(m_total));
left5  = zeros(size(m_total));

for i = 1:601
   
    left1(: , i)  = -gamma * cross(m_total(: , i) , ([B_ext ; 0 ; 0] + [0 ; 0 ; -mu_0 * M_s] .* m_total(: , i)));
    left2(: , i)  = alpha * cross(m_total(: , i) , dmdt(: , i));
    left34(: , i) = -gamma * B_OF * cos(2 * pi * f * t(i)) * cross([sin(phi_B) ; cos(phi_B) ; 0] , m_total(: , i));
    left5(: , i)  = -gamma * B_DL * cos(2 * pi * f * t(i)) * cross(m_total(: , i) , cross([sin(phi_B) ; cos(phi_B) ; 0] , m_total(: , i)));
    
end

-(gamma*cos(phi_B)*(B_DL*omega*1i + B_OF*B_ext*gamma + B_OF*alpha*omega*1i + B_OF*M_s*gamma*mu_0))/...
    (B_ext^2*gamma^2 + B_ext*alpha*gamma*omega*2i + M_s*mu_0*B_ext*gamma^2 - alpha^2*omega^2 + M_s*mu_0*alpha*gamma*omega*1i - omega^2)
-(gamma*cos(phi_B)*(B_DL*B_ext*gamma - B_OF*omega*1i + B_DL*alpha*omega*1i))/...
    (B_ext^2*gamma^2 + B_ext*alpha*gamma*omega*2i + M_s*mu_0*B_ext*gamma^2 - alpha^2*omega^2 + M_s*mu_0*alpha*gamma*omega*1i - omega^2)