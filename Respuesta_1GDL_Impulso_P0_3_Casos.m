% Análisis Dinámico de Estructuras
% Prof. John Esteban Ardila González
% Respuesta de un sistema de 1GDL para un pulso definido por P(t)
clc, clear all, close all

%% Datos de entrada:
k = 1000; % rígidez, N/m
m = 100; % masa, kg
w = (k/m)^0.5; % frecuencia angular, rad/s
T = 2*pi/w; % periodo, s
P0 = 2000; % carga inicial y ctte en t, N
qst0 = P0/k; % desplazamiento estático inicial, m

%% Respuesta Caso 1:q = Po/k [1-cos(wt)] = qsto [1-cos(wt)]:
% Carga Po es ctte desde 0 hasta infinito
t0 = 0; tf = 15; dt = 0.001; % s
t = (t0:dt:tf)'; % vector de tiempo, s
unos = ones(length(t),1); % vector de unos
q = qst0*(unos-cos(w*t)); % respuesta del sistema, desplazamiento, m
qmxab = max(abs(q)); % este es el desplazamiento máximo, m
Aq = qmxab/qst0; % qué tanto se amplifica el desplazamiento estático

%% Respuesta Caso 2:q = qsto[cos(w(t-td))-cos(wt)]
% Po interrumpida en td
td = 3.83*T; % es el tiempo en el que se interrumpe la carga Po, s
disp('Caso 2')
disp(['T/t_d = ',num2str(T/td)])
disp(['t_d/t_f = ',num2str(td/tf)])

for i=1:length(t)
    if t(i)<td
        q2(i) = q(i); % respuesta en desplazamiento para Caso 1
    else
        q2(i) = qst0*(cos(w*(t(i)-td))-cos(w*t(i))); % respuesta para Caso 2 después de td
    end
end
qmxab2 = max(abs(q2));
Aq2 = qmxab2/qst0;

%% Respuesta Caso 3:q = qst0(1-cos(wt))+qst0/td(sen(wt)/w+t)
td = 4*T; % es el tiempo en el que se interrumpe la carga Po, s
disp('Caso 3')
disp(['T/t_d = ',num2str(T/td)])
disp(['t_d/t_f = ',num2str(td/tf)])

qd = qst0*(1-cos(w*td))+qst0/td*(sin(w*td)/w-td);
dqd = qst0*(w*sin(w*td)+1/td*(cos(w*td)-1));

for i=1:length(t)
    if t(i)<td
        q3(i) = qst0*(1-cos(w*t(i)))+qst0/td*(sin(w*t(i))/w-t(i)); % respuesta en desplazamiento para Caso 3 antes de td
    else
        q3(i) = qd*cos(w*t(i))+dqd/w*sin(w*t(i)); % respuesta para Caso 3 después de td
    end
end

%% Gráfica:
figure
plot(t,q,'-r',t,q2,'--k',t,q3,':b','LineWidth',1.4)
xlabel('t (s)'), ylabel('Desplazamiento (m)')
grid on
legend('Caso 1','Caso2','Caso 3')
