% Gradient Descent with Constant Step Size
% Verifying that the set of linear equalities hold for the optimal Q,
% sigmas and t

syms L mu gamma rho
alpha = L/(2*(L-mu));

% optimal Q, t, sigmas
Q = zeros(8,8);
Q = sym(Q);
Q(5,5) = rho^2 + (2*L*gamma*mu*rho)/(L - mu) - 1;
Q(5,6) = -rho^2 - (2*L*gamma*mu*rho)/(L - mu) + 1;
Q(5,7) = (gamma*(mu - L + L*rho + mu*rho))/(L - mu);
Q(6,6) = rho^2 + (2*L*gamma*mu*rho)/(L - mu) - 1;
Q(6,7) = -(gamma*(mu - L + L*rho + mu*rho))/(L - mu);
Q(7,7) = - gamma^2 + (2*rho*gamma)/(L - mu);
for i = 5:8
    for j = i+1:8
        Q(j,i) = Q(i,j);
    end
end

t = rho^2;
sig2 = 2*gamma*rho;
sig5 = 2*gamma*rho;

syms fs fk fp xs xk gk gp
zbar = [1; fs; fk; fp; xs; xk; gk; gp];
zbarT = [1 fs fk fp xs xk gk gp];

% constraints (h >= 0)
h2 = fk - fs - alpha*(1/L*gk^2 + mu*(xk-xs)^2 + 2*mu/L*gk*(xs-xk));
h5 = fs - fk - gk*(xs - xk) - alpha*(1/L*gk^2 + mu*(xs-xk)^2 - 2*mu/L*gk*(xk-xs));

% check whether linear equalities hold
pz = t*(xk - xs)^2 - (xk - gamma*gk - xs)^2;
check = simplify(zbarT*Q*zbar + sig2*h2 + sig5*h5);
diff = simplify(pz - check); % should be zero