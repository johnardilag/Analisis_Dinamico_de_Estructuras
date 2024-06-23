% Análisis Dinámico de Estructuras
% Prof. John Esteban Ardila González
% Ejercicio de aplicación #2
clc, clear all, close all

%% Valores de entrada:
m = 100; % masa en kg
k = 1000; % rigidez en N/m
w = (k/m)^0.5; % fecuencia angular en rad/s
ccr = 2*m*w; % coeficiente de amortiguamiento crítico kg/s
zeta = 5/100; % razón de amortiguamiento
c = zeta*ccr; % coeficiente de amortiguamiento en kg/s
wD = w*(1-zeta^2)^0.5; % frecuencia amortiguada en rad/s
Tn = 2*pi/w; % periodo en s
TD = 2*pi/wD; % periodo amortiguado en s

% Condiciones iniciales:
q0 = 0.1; % desplazamiento inicial en m
dq0 = 1; % velocidad inicial en m/s

%% Respuesta del sistema:
t = (0:0.01:10)'; % vector de tiempo desde 0 s hasta 5 s @ 0.01 s
uno = ones(length(t),1); % vector de unos
qmax = (((dq0+q0*zeta*w)/(wD))^2+(q0)^2)^0.5;
theta = atan((dq0+q0*zeta*w)/(wD*q0));
q = qmax*cos(wD*t-uno*theta).*exp(-zeta*w*t); % desplazamiento en función del tiempo
dq = -qmax*(wD*sin(wD*t-uno*theta)+zeta*w*cos(wD*t-uno*theta)).*exp(-zeta*w*t); % velocidad en función del tiempo
ddq = qmax*(zeta^2*w^2+wD^2)*cos(wD*t-uno*theta).*exp(-zeta*w*t); % aceleración en función del tiempo

%% Energía del sistema:
ET = 1/2*m*dq.^2; % energía cinética
EV = 1/2*k*q.^2; % energía potencial
EM = (ET + EV); % energía mecánica
ED = (uno*max(EM) - (EM)); % energía disipada


%% Gráficas:
% Desplazamiento, velocidad y aceleración
figure
subplot(311), plot(t,q,'-k')
grid on, xlabel('t (s)'), ylabel('q (m)')
subplot(312), plot(t,dq,'-k')
grid on, xlabel('t (s)'), ylabel('dq/dt (m/s)')
subplot(313), plot(t,ddq,'-k')
grid on, xlabel('t (s)'), ylabel('d^2q/dt^2 (m/s^2)')

%Energía cinética, potencial y mecánica
figure
plot(t,ET,'-r',t,EV,'-b',t,EM,'-k',t,ED,'-m')
legend('Energía Cinética','Energía Potencial','Energía Mecánica','Energía Disipada')
grid on, xlabel('t (s)'), ylabel('Energía (J)')

% Velocidad/frecuencia angular VS Desplazamiento
figure
comet(q,dq/w)
xlabel('q (cm)'), ylabel('(dq/dt)/\omega (m)')
grid on

disp(['w = ',num2str(w),' rad/s'])
disp(['w_D = ',num2str(wD),' rad/s'])
disp(['T = ',num2str(Tn),' s'])
disp(['T_D = ',num2str(TD),' s'])
disp(['q_max = ',num2str(qmax),' m'])
