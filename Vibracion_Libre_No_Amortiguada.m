% Análisis Dinámico de Estructuras
% Prof. John Esteban Ardila González
% Ejercicio de aplicación #1
clc, clear all, close all

%% Valores de entrada:
m = 100; % masa en kg
k = 1000; % rigidez en N/m
w = (k/m)^0.5; % fecuencia angular en rad/s
T = 2*pi/w; % periodo en s

% Condiciones iniciales:
q0 = 0.1; % desplazamiento inicial en m
dq0 = 0.0; % velocidad inicial en m/s

%% Respuesta del sistema:
t = (0:0.01:6)'; % vector de tiempo desde 0 s hasta 5 s @ 0.01 s
q = q0*cos(w*t) + dq0/w*sin(w*t); % desplazamiento en función del tiempo
dq = -q0*w*sin(w*t) + dq0*cos(w*t); % velocidad en función del tiempo
ddq = -q0*w^2*cos(w*t) + dq0*w*sin(w*t); % aceleración en función del tiempo

%% Energía del sistema:
ET = 1/2*m*dq.^2; % energía cinética
EV = 1/2*k*q.^2; % energía potencial
EM = ET + EV; % energía total
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
plot(t,ET,'-r',t,EV,'-b',t,EM,'-k')
legend('Energía Cinética','Energía Potencial','Energía Mecánica')
grid on, xlabel('t (s)'), ylabel('Energía (J)')

% Velocidad/frecuencia angular VS Desplazamiento
figure
comet(q,dq/w)
xlabel('q (cm)'), ylabel('(dq/dt)/\omega (m)')
grid on