% Análisis Dinámico de Estructuras
% Prof. John Esteban Ardila González
% Respuesta de un sistema lineal de 1-GDL sometido a un sísmo
% Método de Newmark (pag. 168, Chopra, 5th Edition)
clc, clear all, close all

%% Datos de entrada
nsis = 'Loma_Prieta_1989';
g = 9.81; % aceleración de la gravedad en m/s^2
m = 100/g; % masa del sistema de 1-GDL en kg
k = 100; % rigidez del sistema de 1-GDL en N/m
zeta = 5/100; % coeficiente de amortiguamiento en %
DS = load([nsis,'.dat']); % datos del sismo [t(s) ddug (g)]
ddug = DS(:,2); % aceleración del suelo como fracción de la g
t = DS(:,1); % vector de tiempo en s
nD = length(t); % tamaño del vector de tiempo
dt = t(2); % paso del tiempo en s, depende del registro
MN = 1;  % Método de Newmark: (1) Media o (2) Lineal
c = 2*m*(k/m)^0.5*zeta; % coeficiente de amortiguamiento en kg/s o N.s/m
T = 2*pi*(m/k)^0.5; % Periodo fundamental de vibración, en s

%% Aplicación del Método de Newmark
if MN == 1
    gamma = 1/2; beta = 1/4; % método de aceleración media
else
    gamma = 1/2; beta = 1/6; % método de aceleración lineal
end

% Cálculos iniciales
a1 = 1/(beta*dt^2)*m + gamma/(beta*dt)*c;
a2 = 1/(beta*dt)*m + (gamma/beta-1)*c;
a3 = (1/(2*beta)-1)*m + dt*(gamma/(2*beta)-1)*c;
K = k + a1;
% Preasignación: Condición inicial para q(1) = dq(1) = ddq(1) = 0 
q = zeros(nD,1); dq = zeros(nD,1); ddq = zeros(nD,1);
ddug(end+1) = 0; % extensión del vector hasta n+1

for i=1:nD-1
     p = -m*ddug(i+1)*g + a1*q(i) + a2*dq(i) + a3*ddq(i);
     q(i+1) = p/K;
     dq(i+1) = gamma/(beta*dt)*(q(i+1)-q(i)) + (1-gamma/beta)*dq(i) + dt*(1-gamma/(2*beta))*ddq(i);
     ddq(i+1) = 1/(beta*dt^2)*(q(i+1)-q(i)) - 1/(beta*dt)*dq(i) - (1/(2*beta)-1)*ddq(i);
end
ddug(end) = []; % reducción del vector hasta n

% Valores máximos absolutos de desplazamiento, velocidad y aceleración (respuesta y sismo):
[qmax, tqmax] = max(abs(q)); % Buscar desplazamiento máximo absoluto
[dqmax,tdqmax] = max(abs(dq)); % Buscar velocidad máxima absoluta
[ddqmax,tddqmax] = max(abs(ddq)); % Buscar aceleración máxima absoluta
[ddumax,tddumax] = max(abs(ddug)); % Buscar aceleración máxima absoluta (sismo)

%% Gráficas de la respuesta
tdesf = 1; % desfase para el texto en t
lw = 1.3; % ancho de línea

% Gráfica del registro sísmico
figure
plot(t,ddug,'LineWidth',lw)
hold on
plot(t(tddumax),ddug(tddumax),'or','LineWidth',lw)
text(t(tddumax)+tdesf,ddug(tddumax),[num2str(ddug(tddumax),'%.2f'),' g'])
hold off
xlabel('t (s)'), ylabel('ddu_g (g)')
grid on
title(replace(nsis,'_',' '))

% Gráficas de respuesta del sistema de 1GDL
figure
sgtitle(['T = ',num2str(T,'%.2f'),' s; \zeta = ',num2str(zeta*100),'%'])
subplot(311)
plot(t,q,'LineWidth',lw)
hold on
plot(t(tqmax),q(tqmax),'or','LineWidth',lw)
text(t(tqmax)+tdesf,q(tqmax),[num2str(q(tqmax),'%.2f'),' m'])
hold off
xlabel('tiempo (s)'), ylabel('Desplazamiento (m)')
grid on
title('Histórico de Desplazamiento')

subplot(312)
plot(t,dq,'LineWidth',lw)
hold on
plot(t(tdqmax),dq(tdqmax),'or','LineWidth',lw)
text(t(tdqmax)+tdesf,dq(tdqmax),[num2str(dq(tdqmax),'%.2f'),' m/s'])
hold off
xlabel('tiempo (s)'), ylabel('Velocidad (m/s)')
grid on
title('Histórico de Velocidad')

subplot(313)
plot(t,ddq/g,'LineWidth',lw)
hold on
plot(t(tddqmax),ddq(tddqmax)/g,'or','LineWidth',1.4)
text(t(tddqmax)+tdesf,ddq(tddqmax)/g,[num2str(ddq(tddqmax)/g,'%.2f'),' g'])
hold off
xlabel('tiempo (s)'), ylabel('Aceleración (g)')
grid on
title('Histórico de Aceleración')
