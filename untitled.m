% Definindo a matriz A (exemplo arbitrário)
A = [-1 2; -3 -4];

% Definindo a variável da LMI: P é uma matriz simétrica
n = size(A,1); % Tamanho da matriz A
P = sdpvar(n,n,'symmetric'); % Variável simétrica no YALMIP

% Definindo as restrições da LMI
LMI1 = P >= eye(n); % P deve ser definida positiva (P - I >= 0)
LMI2 = A'*P + P*A <= -eye(n); % Condição de estabilidade (negativa definida)

% Combinando as restrições
Constraints = [LMI1, LMI2];

% Resolvendo o problema com o MOSEK
optimize(Constraints, [], sdpsettings('solver','sdpt3'));

% Extraindo o resultado
P_value = value(P);
disp('Matriz P encontrada:');
disp(P_value);