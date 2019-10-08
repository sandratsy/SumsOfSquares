% Gradient Descent with Wolfe Conditions
% Verifying that equation (61) holds for the optimal sigmas and t
syms mu L const1 const2
alpha = L/(2*(L-mu));
syms fs fk fp xs xk xp gk gp

% optimal t, sigmas and thetas
t = 2*mu*const1*(1-const2)/L;
sigma5 = 2*mu*const1*(1-const2)/L;
sigma7 = 1;

% constraints (h >= 0)
h5 = fs - fk - gk*(xs-xk) - alpha*(1/L*gk^2 + mu*(xs-xk)^2 - 2*mu/L*gk*(xk-xs));
h7 = fk - fp - const1*(1-const2)/L*gk^2;

% check whether equation (61) holds
pz = (1 - t)*(fk - fs) - (fp - fs);
sosterm = const1*(1-const2)/(L-mu)*(gk + mu*(xs-xk))^2;
diff = simplify(pz - (sosterm + sigma5*h5 + sigma7*h7)); % should be zero