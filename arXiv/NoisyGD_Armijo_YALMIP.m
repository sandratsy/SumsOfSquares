% (Noisy) Gradient Descent with Armijo Rule
% Searching for a degree-1 certificate using YALMIP
mu = 3; L = 5;
alpha = L/(2*(L-mu));
delta = 0.5;
eta = 1.01; eps = 0.2; % requirement: eps < (1-delta)/(1+delta)^2;

sdpvar fs fk fp xs xk xp gk gp t
obj = t*(fk - fs) - (fp - fs);

% inequality constraints (c >= 0)
c1 = fk - fp - gp*(xk-xp) - alpha*(1/L*(gk-gp)^2 + mu*(xk-xp)^2 - 2*mu/L*(gp-gk)*(xp-xk));
c2 = fk - fs - alpha*(1/L*gk^2 + mu*(xk-xs)^2 + 2*mu/L*gk*(xs-xk));
c3 = fp - fk - gk*(xp-xk) - alpha*(1/L*(gp-gk)^2 + mu*(xp-xk)^2 - 2*mu/L*(gk-gp)*(xk-xp));
c4 = fp - fs - alpha*(1/L*gp^2 + mu*(xp-xs)^2 + 2*mu/L*gp*(xs-xp));
c5 = fs - fk - gk*(xs-xk) - alpha*(1/L*gk^2 + mu*(xs-xk)^2 - 2*mu/L*gk*(xk-xs));
c6 = fs - fp - gp*(xs-xp) - alpha*(1/L*gp^2 + mu*(xs-xp)^2 - 2*mu/L*gp*(xp-xs));
c7 = fk - fp - (2*eps)/(eta*L)*((1-delta)/(1+delta)^2 - eps)*(1-delta)^2*gk^2;

% define polynomials
[s1, v1] = polynomial([fs fk fp xs xk xp gk gp],0);
[s2, v2] = polynomial([fs fk fp xs xk xp gk gp],0);
[s3, v3] = polynomial([fs fk fp xs xk xp gk gp],0);
[s4, v4] = polynomial([fs fk fp xs xk xp gk gp],0);
[s5, v5] = polynomial([fs fk fp xs xk xp gk gp],0);
[s6, v6] = polynomial([fs fk fp xs xk xp gk gp],0);
[s7, v7] = polynomial([fs fk fp xs xk xp gk gp],0);

% feasible set
F = [sos(obj-[s1 s2 s3 s4 s5 s6 s7]*[c1; c2; c3; c4; c5; c6; c7]), ...
    sos(s1), sos(s2),sos(s3),sos(s4),sos(s5),sos(s6),sos(s7)];
[sol,v,Q] = solvesos(F,t,[],[v1;v2;v3;v4;v5;v6;v7;t]);