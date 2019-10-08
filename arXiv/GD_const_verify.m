% Gradient Descent with Constant Step Size
% Verifying that equation (65) holds for the optimal sigmas and t
syms L mu gamma rho
alpha = L/(2*(L-mu));
syms fs fk fp xs xk gk gp

% optimal t, sigmas
t = rho^2;
sig2 = 2*gamma*rho;
sig5 = 2*gamma*rho;

% constraints (h >= 0)
h2 = fk - fs - alpha*(1/L*gk^2 + mu*(xk-xs)^2 + 2*mu/L*gk*(xs-xk));
h5 = fs - fk - gk*(xs - xk) - alpha*(1/L*gk^2 + mu*(xs-xk)^2 - 2*mu/L*gk*(xk-xs));

pz = t*(xk - xs)^2 - (xk - gamma*gk - xs)^2;

% check whether (65) holds in 1st regime: when gamma in (0, 2/(L+mu)]
assume(rho == 1-gamma*mu);
sosterm1 = gamma*(2-gamma*(L+mu))/(L-mu)*(gk-mu*(xk-xs))^2;
diff1 = simplify(pz - (sosterm1 + sig2*h2 + sig5*h5)); % should be zero

% check whether (65) holds in 2nd regime: when gamma in [2/(L+mu), 2/L)
assume(rho == gamma*L - 1);
sosterm2 = gamma*(gamma*(L+mu)-2)/(L-mu)*(gk-L*(xk-xs))^2;
diff2 = simplify(pz - (sosterm2 + sig2*h2 + sig5*h5)); % should be zero