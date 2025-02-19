clc, clear;
% Iniciar YALMIP
yalmip('clear')

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
A_p = [0, 1;0, 0];
B_p = [0, 0; m1^-1, -m2^-1];
C_p = [1, 0];
D_p = [0, 0];

%% Sistema de malha fechada e formulação de problemas
% Matrizes do controlador

A_c = [-1.7321, 1;-1.0014, -0.0532];
B_c = [1.7321; 1];
C_c = [-0.7071, -26.6009; 0.7071, 26.6009];
D_c = [0; 0];

Bw = [0; m1^-1];

Mn1 = [1, -1, -1, 1];
Mn2 = [1, -1, -1, 1];
% Matriz de influência
Mn = [Mn1, zeros(1,4); zeros(1,4), Mn2];

N1 = [1, 1, -1; 1, 0, 0; 0, 1, 0; 0, 0, 1];

N2 = [1, 1, -1; 1, 0, 0; 0, 1, 0; 0, 0, 1];

N = [N1, zeros(4,3); zeros(4,3), N2];

% Matriz pseudo-inversa

Mps_inv = 0.25 * [Mn1', zeros(4,1); zeros(4,1), Mn2'];

W = [100, 0, 0, 0, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0, 0, 0; 
     0, 0, 1, 0, 0, 0, 0, 0; 
     0, 0, 0, 1, 0, 0, 0, 0;
     0, 0, 0, 0, 1, 0, 0, 0;
     0, 0, 0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 0, 0, 1];

A0 = [A_p + B_p * D_c * C_p, B_p  * C_c; B_c * C_p, A_c];

A = [A0, zeros(4,6); zeros(6,4), zeros(6,6)];
B_ = [B_p; zeros(8,2)];
M = Mn;
B = B_ * M;

%C = [Mps_inv * D_c * C_p, Mps_inv * C_c, N];
C = [Mps_inv * D_c * C_p, Mps_inv * C_c, zeros(m_a, n_f)]; % Garanta 8x10

C_ = N' * W * C;
Bw_ = [Bw; 0; 0];
L_f = [zeros(n_p,n_f); zeros(n_c,n_f); eye(n_f)];

L_c = [zeros(n_p,n_c); eye(n_c); zeros(n_f,n_c)];

%L = [ L_c, L_f];
L = [zeros(n_p, m_a); eye(m_a); zeros(n_f, m_a)]; % n = 10, m_a = 8


P_ = sdpvar(n, n, 'symmetric');
J0 = sdpvar(n_p + n_c, n_p + n_c);
Jf = sdpvar(n_f, n);
Kf_ = sdpvar(n_f, n_f);
Ke = sdpvar(n_c + n_f, m_a);
G_ = sdpvar(m_a, n);
S = sdpvar(m_a, m_a, 'diagonal');

P0 = sdpvar(n, n, 'symmetric');

J_ = [null(C_) * J0', Jf'];

%J_ = sdpvar(n, n, 'symmetric');

Z = [zeros(n_p + n_c,n_p + n_c), zeros(n_p + n_c,n_f); zeros(n_f,n_p + n_c), Kf_];

lambda = sdpvar(1, 1);
gamma = sdpvar(1, 1);
mu = sdpvar(1, 1);

delta = 1e-6;

I = eye(8);

Psi_12 = P_ + A * J_' + Z - J_;
He = A * J_' + Z;
Psi_22 = 0.5 * (He + He');

Psi_13 = B * S + L * Ke;
Psi_23 = Psi_13 - G_' - J_ * C';

Psi_a = [-J_ - J_', Psi_12; Psi_12', Psi_22];
Psi_b = [Psi_13, zeros(n,m_a); Psi_23, J_ * C' * W^0.5];
Psi_c = [-2 * S, S * W^0.5; (S * W^0.5)', -gamma * I];


Psi = [Psi_a, Psi_b; Psi_b', Psi_c];

% Psi_w = [Psi, [Bw_; Bw_; zeros(14,1); zeros(14,1)];
%         [Bw_', Bw_', zeros(1,14), zeros(1,14)], -R];

Bw_extended = [Bw_; zeros(n - n_p - n_c, 1)]; % n = 10
Psi_w = [Psi, Bw_extended; 
         Bw_extended', -R];

% Definir u_bar como um vetor 8x1 (50 mN para cada atuador)

u_ = 50e-3 * ones(m_a, 1);

% Restrições


% Função objetivo
objective = rho1 * lambda + rho2 * gamma + rho3 * mu;

% Solver Configurações
options = sdpsettings('solver', 'sedumi', 'verbose', 2, 'debug', 1, 'cachesolvers', 1);

% Otimização
sol = optimize(constraints, objective, options);

% Resultados
% Restrições Corrigidas
constraints = [];
constraints = [constraints, Psi_w <= -delta * eye(size(Psi_w))];
constraints = [constraints, P0 <= lambda * eye(size(P0))];
constraints = [constraints, sigma - mu >= delta];
constraints = [constraints, [J_ + J_' - P_, eye(n); eye(n), P0] >= delta * eye(2*n)]; % Lyapunov

for i = 1:m_a
    constraints = [constraints, [mu * u_(i)^2, G_(i,:); G_(i,:)', P_] >= delta * eye(n+1)];
end

if sol.problem == 0
    disp('Problema resolvido com sucesso.');
    disp('Solução ótima:');
    disp(value(objective));
    disp('Variáveis de decisão:');
    disp('P_:'); disp(value(P_));
else
    disp('Erro na solução do problema:');
    disp(sol.info);
    disp('Verifique as restrições ou inicializações.');
end