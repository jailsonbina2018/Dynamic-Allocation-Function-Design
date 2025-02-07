% Simulando o sistema da imagem em MATLAB

% Dimensões dos sistemas
n_p = 2;
m_c = 2;
n_c = 2;
m_a = 8;
n_w = 2;

n_f = m_a - n_c;
n = n_p + n_c + n_f;

% Parâmetros conhecidos
m1 = 1000;
m2 = 1000;
R = 1;
sigma = 1;
rho1 = 2;
rho2 = 0.15;
rho3 = 1000;

% Matrizes da planta
A_p = [0, 1; 0, 0];
B_p = [0, 0; m1^-1, -m2^-1];
C_p = [1, 0];
D_p = [0, 0];

%% Sistema de malha fechada e formulação de problemas
% Matrizes do controlador
A_c = [-1.7321, 1; -1.0014, -0.0532];
B_c = [1.7321; 1];
C_c = [-0.7071, -26.6009; 0.7071, 26.6009];
D_c = [0; 0];

Bw = [0; m1^-1];

Mn1 = [1, -1, -1, 1];
Mn2 = [1, -1, -1, 1];
% Matriz de influência
Mn = [Mn1, zeros(1, 4); zeros(1, 4), Mn2];

N1 = [1, 1, -1; 1, 0, 0; 0, 1, 0; 0, 0, 1];
N2 = [1, 1, -1; 1, 0, 0; 0, 1, 0; 0, 0, 1];
N = [N1, zeros(4, 3); zeros(4, 3), N2];

% Matriz pseudo-inversa
Mps_inv = 0.25 * [Mn1', zeros(4, 1); zeros(4, 1), Mn2'];

W = [100, 0, 0, 0, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0, 0, 0;
     0, 0, 0, 1, 0, 0, 0, 0;
     0, 0, 0, 0, 1, 0, 0, 0;
     0, 0, 0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 0, 0, 1];

A0 = [A_p + B_p * D_c * C_p, B_p * C_c; B_c * C_p, A_c];
A = [A0, zeros(4, 6); zeros(6, 4), zeros(6, 6)];
B_ = [B_p; zeros(8, 2)];
M = Mn;
B = B_ * M;

C = [Mps_inv * D_c * C_p, Mps_inv * C_c, N];
C_ = N' * W * C;
Bw_ = [Bw; zeros(4, 1); zeros(4, 1)];
L_f = [zeros(n_p, n_f); zeros(n_c, n_f); eye(n_f)];
L_c = [zeros(n_p, n_c); eye(n_c); zeros(n_f, n_c)];
L = [L_c, L_f];

E_c = [-0.1708, -0.0059, -0.0032, -0.0053, 0.0111, -0.0066, -0.0066, 0.0066;
       -0.0895, -0.0118, -0.0103, 0.0057, -0.0024, 0.0048, 0.0048, -0.0048];

E_f = [-0.0480, -0.0103, 0.1031, -0.0662, 0.0747, 0.1042, 0.1042, -0.1042;
       -0.0805, 0.0057, -0.0661, 0.0225, -0.0512, -0.0537, -0.0537, 0.0537;
        0.1277, -0.0018, 0.0566, -0.0389, -0.0080, 0.0722, 0.0722, -0.0722;
       -0.0757, 0.0037, 0.0790, -0.0408, 0.0722, 0.9615, -0.3662, 0.3662;
       -0.0757, 0.0037, 0.0790, -0.0408, 0.0722, -0.3662, 0.9615, 0.3662;
        0.0757, -0.0037, -0.0790, 0.0408, -0.0722, 0.3662, 0.3662, 0.9615];

K_f = [-1.8960, 0.9476, -0.9437, 0.0053, 0.0053, -0.0053;
        0.9103, -1.8750, -0.9090, 0.0013, 0.0013, -0.0013;
       -0.8741, -0.8266, -1.8904, -0.0007, -0.0007, 0.0007;
       -0.0020, -0.0009, 0.0118, -1.6716, 0.5010, -0.5010;
       -0.0020, -0.0009, 0.0118, 0.5010, -1.6716, -0.5010;
        0.0020, 0.0009, -0.0118, 0.5010, -1.6716, -0.5010];

E = [E_c', E_f'];

x_p = [-0.18; 0];
x_c = [0; 0];
x_f = zeros(n_f, 1);

x = [x_p; x_c; x_f];
  
% Tempo de simulação
dt = 0.01; % Passo de tempo
T = 180;   % Duração da simulação
t = 0:dt:T;

% Sinal de entrada w(t)
w = zeros(1, length(t));
w(t >= 0 & t <= 36) = 0.16667;

% Função de saturação e phi
% Definindo os limites de saturação simétricos
u_bar = 0.5; % Limite superior simétrico (50 mN)

% Função de saturação baseada em sign
sat = @(x) sign(x) .* min(abs(x), u_bar);

% Atualizando a função phi
phi = @(yf) sat(yf) - yf;

% Inicialização dos estados
x_states = zeros(n, length(t));
y_outputs = zeros(size(C, 1), length(t)); % Armazena as saídas
x_states(:, 1) = x;

% Simulação do sistema
for k = 1:length(t)-1
    y_f = C * x_states(:, k);
    dx = (A + L_f * K_f * C_) * x_states(:, k) + (B + L * E) * phi(y_f) + Bw_ * w(k);
    x_states(:, k + 1) = x_states(:, k) + dx * dt;
    y_outputs(:, k) = y_f; % Armazena a saída
end

% Plotagem das saídas do sistema
figure;
plot(t, y_outputs');
xlabel('Tempo (s)');
ylabel('Saídas do sistema (y)');
title('Resposta temporal das saídas do sistema');
legend(arrayfun(@(i) sprintf('y%d', i), 1:size(C, 1), 'UniformOutput', false));
grid on;
