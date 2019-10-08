% Gradient Descent with Exact Line Search
% Verifying that equation (42) holds for the optimal sigmas, thetas and t

syms L mu
alpha = L/(2*(L-mu));
kappa = mu/L;
syms fs fk fp xs xk xp gk gp

% optimal t, sigmas and thetas
t = ((L-mu)/(L+mu))^2;
sig1 = (L-mu)/(L+mu);
sig5 = 2*mu*(L-mu)/(L+mu)^2;
sig6 = 2*mu/(L+mu);
the1 = -1;
the2 = -2/(L+mu);

% sos term
q1 = -(1+mu^(1/2)/L^(1/2))^2/(1+kappa)*(xk-xs-gk/(L^(1/2)*mu^(1/2))) + (xp-xs+gp/(L^(1/2)*mu^(1/2)));
q2 = (1-mu^(1/2)/L^(1/2))^2/(1+kappa)*(xk-xs+gk/(L^(1/2)*mu^(1/2))) - (xp-xs-gp/(L^(1/2)*mu^(1/2)));
sos = mu/4*(q1^2/(1+mu^(1/2)/L^(1/2)) + q2^2/(1-mu^(1/2)/L^(1/2)));

% constraints (h >= 0, v == 0)
h1 = fk - fp - gp*(xk - xp) - alpha*(1/L*(gk-gp)^2 + mu*(xk-xp)^2 - 2*mu/L*(gp-gk)*(xp-xk));
h5 = fs - fk - gk*(xs - xk) - alpha*(1/L*gk^2 + mu*(xs-xk)^2 - 2*mu/L*gk*(xk-xs));
h6 = fs - fp - gp*(xs - xp) - alpha*(1/L*gp^2 + mu*(xs-xp)^2 - 2*mu/L*gp*(xp-xs));
v1 = gp*(xp - xk);
v2 = gp*gk;

% check whether equation (42) holds
pz = t*(fk - fs) - (fp - fs);
check = simplify(sos + sig1*h1 + sig5*h5 + sig6*h6 + the1*v1 + the2*v2);
diff = simplify(pz - check); % should be zero