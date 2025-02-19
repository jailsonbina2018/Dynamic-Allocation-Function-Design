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

C = [Mps_inv * D_c * C_p, Mps_inv * C_c, N];
C_ = N' * W * C;
Bw_ = [Bw; 0; 0];
L_f = [zeros(n_p,n_f); zeros(n_c,n_f); eye(n_f)];

L_c = [zeros(n_p,n_c); eye(n_c); zeros(n_f,n_c)];

L = [ L_c, L_f];

P_ = sdpvar(n, n);
P0 = sdpvar(n, n);
J_ = sdpvar(n, n);
S = sdpvar(m_a, m_a, 'symmetric');

G_ = sdpvar(m_a, n);
Kf_ = sdpvar(n_f, n_f);
Ke = sdpvar(m_a, m_a);

Z = [zeros(n_p + n_c,n_p + n_c), zeros(n_p + n_c,n_f); zeros(n_f,n_p + n_c), Kf_];

lambda = sdpvar(1, 1);
gamma = sdpvar(1, 1);
mu = sdpvar(1, 1);

delta = 1e-12;

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

Psi_w = [Psi, [Bw_; Bw_; zeros(14,1); zeros(14,1)];
        [Bw_', Bw_', zeros(1,14), zeros(1,14)], -R];

% Definir u_bar como um vetor 8x1 (50 mN para cada atuador)

u_ = 50e-3 * ones(m_a, 1);

% Restrições
constraints = [];

constraints = [constraints, Psi_w <= delta * eye(size(Psi_w))];
constraints = [constraints, P0 <= lambda * eye(size(P0))];
constraints = [constraints, sigma - mu >= delta];
constraints = [P0, eye(10); eye(10), J_ - J_' - P_] >= delta * eye(size([P0, eye(10); eye(10), J_ - J_' - P_]));

for i = 1:m_a
    constraints = [constraints, [P_, G_(i,:)'; G_(i,:), mu * u_(i,:)^2] >= delta * eye(size([P_, G_(i,:)'; G_(i,:), mu * u_(i,:)^2]))];
end

% Função custo
objective = rho1 * lambda + rho2 * gamma + rho3 * mu;

% Solver
options = sdpsettings('solver', 'sedumi', 'verbose', 2, 'debug', 1);
sol = optimize(constraints, objective, options);

% Resultados
if sol.problem == 0
    disp('Problema resolvido com sucesso.');
    disp('Solução ótima:');
    disp(value(objective));
    disp('Variáveis de decisão:');
    disp('P_:'); disp(value(P_));
    disp('gamma:'); disp(value(gamma));
    disp('mu:'); disp(value(mu));
else
    disp('Erro na solução do problema:');
    disp(sol.info);
end





