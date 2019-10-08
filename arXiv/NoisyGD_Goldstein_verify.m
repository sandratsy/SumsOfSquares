% (Noisy) Gradient Descent with Goldstein Rule
% Verifying that equation (57) holds for the optimal sigmas and t
syms mu L delta eps
alpha = L/(2*(L-mu));
syms fs fk fp xs xk xp gk gp

% optimal t, sigmas and thetas
t = 4*mu*eps*(1-delta)^2/L*((1-delta)/(1+delta)^2 - (1-eps));
sigma5 = 4*mu*eps*(1-delta)^2/L*((1-delta)/(1+delta)^2 - (1-eps));
sigma7 = 1;

% constraints (h >= 0)
h5 = fs - fk - gk*(xs-xk) - alpha*(1/L*gk^2 + mu*(xs-xk)^2 - 2*mu/L*gk*(xk-xs));
h7 = fk - fp - (2*eps*(1-delta)^2)/L*((1-delta)/(1+delta)^2-1+eps)*gk^2;

% check whether equation (57) holds
pz = (1 - t)*(fk - fs) - (fp - fs);
sosterm = 2*eps^2/(L-mu)*(gk + mu*(xs-xk))^2;
diff = simplify(pz - (sosterm + sigma5*h5 + sigma7*h7)); % should be zero