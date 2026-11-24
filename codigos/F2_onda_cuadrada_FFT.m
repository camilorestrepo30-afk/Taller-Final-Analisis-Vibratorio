%% F2 - Excitacion de onda cuadrada y analisis FFT
% m x'' + c x' + k x = F(t) , F(t) = F0*sign(sen(2 pi f t))
clear; clc; close all;

%% ---------- GRAFICAS CON FONDO BLANCO ----------
set(groot, 'defaultFigureColor', 'w');          % fondo de la ventana
set(groot, 'defaultAxesColor', 'w');            % fondo de los ejes
set(groot, 'defaultAxesXColor', 'k', 'defaultAxesYColor', 'k', 'defaultAxesZColor', 'k');
set(groot, 'defaultTextColor', 'k', 'defaultAxesGridColor', [0.15 0.15 0.15]);
try, set(groot, 'defaultLegendColor', 'w', 'defaultLegendTextColor', 'k'); end
set(groot, 'defaultAxesFontSize', 11);
try, set(groot, 'defaultFigureTheme', 'light'); end   % MATLAB R2025a+ (modo oscuro)

%% ---------- PARAMETROS ----------
m  = 1;  k = 355.3;  c = 1.0;      % sistema (fn ~ 3 Hz)
F0 = 10;                          % amplitud de la onda cuadrada [N]
f_list = [1.5 1.0];               % frecuencias a estudiar [Hz] (1.5 = enunciado)
fs = 1000;                        % frecuencia de muestreo [Hz]
t_trans = 40;                     % tiempo para eliminar el transitorio [s]
n_per   = 30;                     % numero de periodos ENTEROS para la FFT

wn = sqrt(k/m);  fn = wn/(2*pi);  zeta = c/(2*sqrt(k*m));
fprintf('fn = %.4f Hz , zeta = %.4f\n', fn, zeta);
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-3);

for f = f_list
    T  = 1/f;
    tf = t_trans + n_per*T;
    t  = 0:1/fs:tf;
    F  = @(t) F0*sign(sin(2*pi*f*t));
    ode = @(t,y) [y(2); (F(t) - c*y(2) - k*y(1))/m];
    [~,y] = ode45(ode, t, [0; 0], opts);
    x = y(:,1)';

    % ---- Estado estacionario: ultimos n_per periodos completos ----
    ie = t > t_trans;  ie(find(ie,1,'last')) = false;   % N = n_per*T*fs muestras
    xs = x(ie);  ts = t(ie);  N = length(xs);

    % ---- FFT unilateral ----
    X  = fft(xs)/N;
    A  = 2*abs(X(1:floor(N/2)+1));  A(1) = A(1)/2;
    fv = (0:floor(N/2))*fs/N;

    % ---- Prediccion teorica por serie de Fourier ----
    % F(t) = (4F0/pi) * sum_{n impar} sen(n w t)/n
    fprintf('\n=== f = %.2f Hz ===\n  n   f_n[Hz]   A_FFT[mm]   A_teo[mm]\n', f);
    for n = 1:7
        Fn = (n/2 ~= round(n/2))*4*F0/(n*pi);       % solo armonicos impares
        r  = n*f/fn;
        At = Fn/k/sqrt((1 - r^2)^2 + (2*zeta*r)^2);
        [~,ib] = min(abs(fv - n*f));
        fprintf('  %d   %6.2f   %9.4f   %9.4f\n', n, n*f, 1e3*A(ib), 1e3*At);
    end
    [Amax, im] = max(A);
    fprintf('  Pico dominante del espectro: %.3f Hz (A = %.3f mm)\n', fv(im), 1e3*Amax);

    figure;
    subplot(2,1,1); plot(ts, 1e3*xs); grid on;
    xlabel('t [s]'); ylabel('x [mm]');
    title(sprintf('Estado estacionario, onda cuadrada f = %.2f Hz', f));
    subplot(2,1,2); stem(fv, 1e3*A, 'filled'); hold on; grid on;
    plot([fn fn], [0 1.1e3*Amax], 'r--');
    xlim([0 10]); xlabel('f [Hz]'); ylabel('|X(f)| [mm]');
    title('Espectro de amplitud (FFT)'); legend('FFT', 'f_n');
end
