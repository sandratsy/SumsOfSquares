% (Noisy) Gradient Descent with Armijo Rule
% Verifying that equation (53) holds for the optimal sigmas and t

syms mu L delta eta eps
alpha = L/(2*(L-mu));
syms fs fk fp xs xk xp gk gp

% optimal t, sigmas and thetas
t = (4*mu*eps)/(eta*L)*((1-delta)/(1+delta)^2 - eps)*(1-delta)^2;
sigma5 = (4*mu*eps)/(eta*L)*((1-delta)/(1+delta)^2 - eps)*(1-delta)^2;
sigma7 = 1;

% constraints (h >= 0)
h5 = fs - fk - gk*(xs-xk) - alpha*(1/L*gk^2 + mu*(xs-xk)^2 - 2*mu/L*gk*(xk-xs));
h7 = fk - fp - (2*eps)/(eta*L)*((1-delta)/(1+delta)^2 - eps)*(1-delta)^2*gk^2;

% check whether equation (53) holds
pz = (1 - t)*(fk - fs) - (fp - fs);
sosterm = (2*eps*(1-delta)^2*(1-delta-eps*(1+delta)^2))/(eta*(L-mu)*(1+delta)^2)*(gk + mu*(xs-xk))^2;
diff = simplify(pz - (sosterm + sigma5*h5 + sigma7*h7)); % should be zero
