% Proximal Gradient Method with constant step size
% Verifying that equation (71) holds for the optimal sigmas and t
syms L mu gamma rho
alpha = L/(2*(L-mu));
syms fs fk fp hs hk hp xs xk gs gk gp sk sp
xp = xk - gamma*(gk + sp);
ss = -gs;

% optimal t and sigmas
t = rho^2;
sig1 = 2*rho/gamma;
sig3 = 2*rho/gamma;
sig7 = 2*rho^2/gamma;
sig9 = 2*rho^2/gamma;

% inequality constraints (h >= 0)
h1 = fk - fp - gp*(xk-xp) - alpha*(1/L*(gk-gp)^2 + mu*(xk-xp)^2 - 2*mu/L*(gp-gk)*(xp-xk));
h3 = fp - fk - gk*(xp-xk) - alpha*(1/L*(gp-gk)^2 + mu*(xp-xk)^2 - 2*mu/L*(gk-gp)*(xk-xp));
h7 = hk - hp - sp*(xk - xp);
h9 = hp - hk - sk*(xp - xk);

pz = t*(gk + sk)^2 - (gp + sp)^2;
% check whether (71) holds in 1st regime: when gamma in (0, 2/(L+mu)]
assume(rho == 1-gamma*mu);
sosterm1 = rho^2*(sk-sp)^2 + (2-gamma*(L+mu))/(gamma*(L-mu))*(gk-gp-mu*gamma*(gk+sp))^2;
diff1 = simplify(pz - (sosterm1 + sig1*h1 + sig3*h3 + sig7*h7 + sig9*h9)); % should be zero

% check whether (71) holds in 2nd regime: when gamma in [2/(L+mu), 2/L)
assume(rho == gamma*L - 1);
sosterm2 = rho^2*(sk-sp)^2 + (gamma*(L+mu)-2)/(gamma*(L-mu))*(gk-gp-L*gamma*(gk+sp))^2;
diff2 = simplify(pz - (sosterm2 + sig1*h1 + sig3*h3 + sig7*h7 + sig9*h9)); % should be zero