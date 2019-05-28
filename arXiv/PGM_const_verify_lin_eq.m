% Proximal Gradient Method with constant step size
% Verifying that the set of linear equalities hold for the optimal Q,
% sigmas and t

syms L mu gamma rho
alpha = L/(2*(L-mu));

% optimal Q, t and lambdas
Q = zeros(14,14);
Q = sym(Q);
Q(11,11) = -(rho*(2*L*gamma + 2*gamma*mu - L*gamma*rho + gamma*mu*rho - 2*L*gamma^2*mu - 2))/(gamma*(L - mu));
Q(11,12) = (rho*(L*gamma + gamma*mu - 2))/(gamma*(L - mu));
Q(11,14) = -(rho*(L + mu - L*rho + mu*rho - 2*L*gamma*mu))/(L - mu);
Q(12,12) = (2*rho)/(gamma*(L - mu)) - 1;
Q(12,14) = rho + (2*mu*rho)/(L - mu) - 1;
Q(13,13) = rho^2;
Q(13,14) = -rho^2;
Q(14,14) = 2*rho^2 + (2*L*gamma*mu*rho)/(L - mu) - 1;
for i = 11:14
    for j = i+1:14
        Q(j,i) = Q(i,j);
    end
end

t = rho^2;
sig1 = 2*rho/gamma;
sig3 = 2*rho/gamma;
sig7 = 2*rho^2/gamma;
sig9 = 2*rho^2/gamma;

syms fs fk fp hs hk hp xs xk gs gk gp sk sp;
zbar = [1; fs; fk; fp; hs; hk; hp; xs; xk; gs; gk; gp; sk; sp];
zbarT = [1 fs fk fp hs hk hp xs xk gs gk gp sk sp];

% Eliminating variables using equality constraints
xp = xk - gamma*(gk + sp);
ss = -gs;

% inequality constraints (h >= 0)
h1 = fk - fp - gp*(xk-xp) - alpha*(1/L*(gk-gp)^2 + mu*(xk-xp)^2 - 2*mu/L*(gp-gk)*(xp-xk));
h3 = fp - fk - gk*(xp-xk) - alpha*(1/L*(gp-gk)^2 + mu*(xp-xk)^2 - 2*mu/L*(gk-gp)*(xk-xp));
h7 = hk - hp - sp*(xk - xp);
h9 = hp - hk - sk*(xp - xk);

% check whether set of linear equalities hold
pz = t*(gk + sk)^2 - (gp + sp)^2;
check = simplify(zbarT*Q*zbar + sig1*h1 + sig3*h3 + sig7*h7 + sig9*h9);
diff = simplify(pz - check); % should be zero