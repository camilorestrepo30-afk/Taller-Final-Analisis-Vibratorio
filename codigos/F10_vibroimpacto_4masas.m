%% F10 - Cadena de 4 masas con topes elasticos (vibro-impacto)
% Diferencias finitas centradas + barrido de frecuencia + FFT (Hann) + THD + Parseval
% pared -k,c- m1 -k,c- m2 -k,c- m3 -k,c- m4 -k,c- pared ;  topes en m2 y m3
clear; clc; close all;

%% ---------- GRAFICAS CON FONDO BLANCO ----------
set(groot, 'defaultFigureColor', 'w');          % fondo de la ventana
set(groot, 'defaultAxesColor', 'w');            % fondo de los ejes
set(groot, 'defaultAxesXColor', 'k', 'defaultAxesYColor', 'k', 'defaultAxesZColor', 'k');
set(groot, 'defaultTextColor', 'k', 'defaultAxesGridColor', [0.15 0.15 0.15]);
try, set(groot, 'defaultLegendColor', 'w', 'defaultLegendTextColor', 'k'); end
set(groot, 'defaultAxesFontSize', 11);
try, set(groot, 'defaultFigureTheme', 'light'); end   % MATLAB R2025a+ (modo oscuro)

%% ================= PARAMETROS (modificables) =================
N     = 4;          % numero de masas
m     = 1;          % masa [kg]
k     = 100;        % rigidez [N/m]
c     = 0.4;        % amortiguamiento [N s/m]
ktope = 5e4;        % rigidez de los topes [N/m]
d     = 0.05;       % holgura de los topes [m]
F0    = 15;         % amplitud de la fuerza sobre m1 [N]
masas_tope = [0 1 1 0]';        % 1 = la masa tiene tope (m2 y m3)

w_v     = linspace(4, 28, 120); % barrido de frecuencia [rad/s]
n_tot   = 200;      % ciclos simulados por frecuencia
n_trans = 100;      % ciclos transitorios descartados
n_est   = 50;       % ultimos ciclos conservados (regimen estacionario)
Nfft    = 2^14;     % muestras para la FFT (potencia de dos)
seg     = 0.5;      % factor de seguridad CFL: dt <= seg * 2/w_max
armonicos = [1 2 3 5];          % multiplos de w a extraer
h_thd     = 2:10;               % armonicos usados en el THD

%% ================= (a) MATRICES M, C, K =================
M = m*eye(N);
K = k*(2*eye(N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1));
C = c*(2*eye(N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1));

% Fuerza de los topes (bilineal con holgura), vectorizada
ftope = @(x) masas_tope .* ( ktope*(x - d).*(x > d) + ktope*(x + d).*(x < -d) );

%% ================= (b) FRECUENCIA MAXIMA Y PASO DE TIEMPO =================
w_lin = sqrt(eig(M\K));                          % sin topes
K_act = K + ktope*diag(masas_tope);              % topes activos (rigidez tangente)
w_act = sqrt(eig(M\K_act));
w_max = max(w_act);
dt_cfl = 2/w_max;
fprintf('Frecuencias naturales SIN topes [rad/s]: %s\n', mat2str(w_lin', 4));
fprintf('Frecuencias naturales CON topes [rad/s]: %s\n', mat2str(w_act', 5));
fprintf('w_max = %.3f rad/s  ->  limite CFL: dt < 2/w_max = %.5f s\n', w_max, dt_cfl);

% Puntos por ciclo = potencia de dos tal que dt = T/ppc <= seg*dt_cfl para TODA w
T_max = 2*pi/min(w_v);
ppc   = 2^nextpow2(T_max/(seg*dt_cfl));
fprintf('Puntos por ciclo ppc = 2^%d = %d\n', log2(ppc), ppc);
fprintf('dt maximo usado (w = %.1f) = %.5f s  ( dt*w_max = %.3f < 2 )\n', ...
        min(w_v), T_max/ppc, T_max/ppc*w_max);
fprintf('Ventana FFT: Nfft = %d muestras = %d ciclos completos\n', Nfft, Nfft/ppc);
if Nfft/ppc > n_est, error('Nfft excede el registro estacionario: reduzca Nfft o aumente n_est'); end

%% ================= RESERVA DE RESULTADOS =================
Nw   = length(w_v);
Aarm = zeros(N, length(armonicos), Nw);   % amplitudes en w, 2w, 3w, 5w
Asub = zeros(N, 2, Nw);                   % amplitudes en w/2 y w/3
A1x2 = zeros(1,Nw);  THD = zeros(1,Nw);
Ek   = zeros(N, Nw);  Ek_t = zeros(N, Nw);
xmax = zeros(N, Nw);  contacto = zeros(1, Nw);
P2map = zeros(6*Nfft/ppc + 1, Nw);        % espectro de x2 hasta 6w (mapa)

whann = 0.5*(1 - cos(2*pi*(0:Nfft-1)/Nfft));   % ventana Hann (periodica)
q     = Nfft/ppc;                               % ciclos en la ventana -> bin de w

tic;
for iw = 1:Nw
    w  = w_v(iw);
    T  = 2*pi/w;
    dt = T/ppc;
    Nt = n_tot*ppc;                  % pasos totales

    %% ---- (c) esquema explicito de diferencias centradas ----
    % M (x+ - 2x + x-)/dt^2 + C (x+ - x-)/(2dt) + K x + ftope(x) = F
    % => x+ = A1 * [ F - K x - ftope(x) + (2M/dt^2) x - B x- ]
    A1 = inv(M/dt^2 + C/(2*dt));
    B  = M/dt^2 - C/(2*dt);
    M2 = 2*M/dt^2;

    Nrec = n_est*ppc;                % muestras del registro estacionario
    i0   = Nt - Nrec + 1;            % primer indice guardado
    Xrec = zeros(N, Nrec + 4);       % +2 muestras a cada lado (derivada 4to orden)

    x_ant = zeros(N,1);  x = zeros(N,1);     % reposo inicial
    for n = 1:(Nt + 2)
        tn = (n-1)*dt;
        F  = [F0*cos(w*tn); zeros(N-1,1)];
        x_sig = A1*(F - K*x - ftope(x) + M2*x - B*x_ant);
        x_ant = x;  x = x_sig;
        j = n - (i0 - 2) + 1;                % posicion dentro de Xrec
        if j >= 1, Xrec(:, j) = x; end
    end
    % x(n) del registro estacionario y velocidad por diferencias de 4to orden
    Xs = Xrec(:, 3:end-2);
    Vs = (-Xrec(:,5:end) + 8*Xrec(:,4:end-1) - 8*Xrec(:,2:end-3) + Xrec(:,1:end-4))/(12*dt);

    % Ventana FFT: ultimas Nfft muestras (numero entero de ciclos)
    Xf = Xs(:, end-Nfft+1:end);
    Vf = Vs(:, end-Nfft+1:end);

    %% ---- (e) FFT con ventana Hann y espectro de potencia unilateral ----
    for i = 1:N
        Xk = fft((Xf(i,:) - mean(Xf(i,:))) .* whann);
        A  = 2*abs(Xk(1:Nfft/2+1))/sum(whann);      % amplitud unilateral
        P  = A.^2/2;                                % potencia unilateral
        for a = 1:length(armonicos)
            b = armonicos(a)*q + 1;
            Aarm(i,a,iw) = max(A(b-1:b+1));
        end
        Asub(i,1,iw) = max(A(round(q/2) + (0:2)));             % w/2
        Asub(i,2,iw) = max(A(floor(q/3) + (1:2)));             % w/3
        if i == 2
            Ah = zeros(size(h_thd));
            for hh = 1:length(h_thd)
                b = h_thd(hh)*q + 1;  Ah(hh) = max(A(b-1:b+1));
            end
            A1x2(iw) = max(A(q:q+2));
            THD(iw)  = 100*sqrt(sum(Ah.^2))/A1x2(iw);
            P2map(:,iw) = P(1:6*q+1)';
        end
    end

    %% ---- (f) energia cinetica por Parseval ----
    % mean(v^2) = (1/N) sum v_n^2 = (1/N^2) sum |V_k|^2
    for i = 1:N
        Vk = fft(Vf(i,:));
        Ek(i,iw)   = 0.5*m*sum(abs(Vk).^2)/Nfft^2;   % dominio de la frecuencia
        Ek_t(i,iw) = 0.5*m*mean(Vf(i,:).^2);         % comprobacion en el tiempo
    end
    xmax(:,iw)   = max(abs(Xs), [], 2);
    contacto(iw) = 100*mean(abs(Xs(2,:)) > d | abs(Xs(3,:)) > d);   % % del tiempo en contacto

    if mod(iw,20) == 0, fprintf('  %3d/%d  w = %.2f rad/s  (%.0f s)\n', iw, Nw, w, toc); end
end
eta = Ek ./ sum(Ek, 1);
fprintf('Error maximo Parseval (frecuencia vs tiempo): %.2e\n', max(abs(Ek(:) - Ek_t(:))./Ek_t(:)));

%% ================= RESULTADOS EN CONSOLA =================
banda = w_v(THD > 5);
if ~isempty(banda)
    fprintf('THD(x2) > 5%% en w de %.2f a %.2f rad/s (posibles sub-bandas, ver grafica)\n', min(banda), max(banda));
end
[thdmax, it] = max(THD);
fprintf('THD maximo de x2 = %.1f %% en w = %.2f rad/s\n', thdmax, w_v(it));
sub = squeeze(Asub(2,:,:)) ./ A1x2;          % fila 1: w/2 , fila 2: w/3
isub = find(max(sub,[],1) > 0.10);           % criterio: > 10% de la fundamental
fprintf('Subarmonicos (>10%% de la fundamental) en %d de %d frecuencias', length(isub), Nw);
if ~isempty(isub), fprintf(', entre %.2f y %.2f rad/s', w_v(isub(1)), w_v(isub(end))); end
fprintf('\n');
[smax, js] = max(sub(:));  [fila, col] = ind2sub(size(sub), js);
nombres = {'w/2', 'w/3'};
fprintf('Subarmonico mas fuerte: %s con %.2f veces la fundamental en w = %.2f rad/s\n', ...
        nombres{fila}, smax, w_v(col));
fprintf('Fraccion media de energia en masas centrales (eta2+eta3): %.3f\n', mean(eta(2,:)+eta(3,:)));

%% ================= GRAFICAS =================
% 1) Amplitud de armonicos de cada masa
figure;
for i = 1:N
    subplot(N,1,i);
    semilogy(w_v, squeeze(Aarm(i,:,:))', '.-'); grid on;
    ylabel(sprintf('A_{x%d} [m]', i));
    if i == 1
        title('Amplitud de la fundamental y superarmonicos');
        legend('\omega', '2\omega', '3\omega', '5\omega');
    end
end
xlabel('\omega [rad/s]');

% 2) THD de x2 y zonas de contacto con los topes
figure;
subplot(2,1,1);
plot(w_v, THD, 'b.-'); grid on; ylabel('THD x_2 [%]');
title('Distorsion armonica total de x_2(t)');
subplot(2,1,2);
plot(w_v, xmax(2,:), 'r.-', w_v, xmax(3,:), 'm.-'); hold on; grid on;
plot(w_v, d*ones(size(w_v)), 'k--');
plot(w_v, contacto/100*max(xmax(:)), 'g-');
xlabel('\omega [rad/s]'); ylabel('max|x_i| [m]');
legend('x_2', 'x_3', 'holgura d', '% tiempo en contacto (escalado)');

% 3) Subarmonicos
figure;
plot(w_v, sub(1,:), 'b.-', w_v, sub(2,:), 'r.-'); grid on;
xlabel('\omega [rad/s]'); ylabel('A / A(\omega)');
legend('\omega/2', '\omega/3'); title('Indicador de subarmonicos en x_2');

% 4) Mapa espectral de x2
figure;
imagesc(w_v, (0:6*q)/q, 10*log10(P2map + 1e-16)); axis xy; colorbar;
xlabel('\omega [rad/s]'); ylabel('f / f_\omega (orden armonico)');
title('Espectro de potencia de x_2 [dB]');

% 5) Fraccion de energia cinetica
figure;
plot(w_v, eta', '.-'); grid on;
xlabel('\omega [rad/s]'); ylabel('\eta_i = E_{k,i}/\Sigma E_{k,j}');
legend('m_1', 'm_2 (tope)', 'm_3 (tope)', 'm_4');
title('Distribucion de energia cinetica (Parseval)');
