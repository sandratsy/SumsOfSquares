% Gradient Descent with Constant Step Size
% Searching for a degree-1 certificate using YALMIP
mu = 3; L = 5;
alpha = L/(2*(L-mu));
gamma = 0.15; % set to any value in (0, 2/L)

sdpvar fs fk fp xs xk gk gp t
xp = xk - gamma*gk; % we substitute x_{k+1} with this equation

% p(z)
obj = t*(xk - xs)^2 - (xp - xs)^2;

% inequality constraints (h >= 0)
h1 = fk - fp - gp*(xk-xp) - alpha*(1/L*(gk-gp)^2 + mu*(xk-xp)^2 - 2*mu/L*(gp-gk)*(xp-xk));
h2 = fk - fs - alpha*(1/L*gk^2 + mu*(xk-xs)^2 + 2*mu/L*gk*(xs-xk));
h3 = fp - fk - gk*(xp-xk) - alpha*(1/L*(gp-gk)^2 + mu*(xp-xk)^2 - 2*mu/L*(gk-gp)*(xk-xp));
h4 = fp - fs - alpha*(1/L*gp^2 + mu*(xp-xs)^2 + 2*mu/L*gp*(xs-xp));
h5 = fs - fk - gk*(xs-xk) - alpha*(1/L*gk^2 + mu*(xs-xk)^2 - 2*mu/L*gk*(xk-xs));
h6 = fs - fp - gp*(xs-xp) - alpha*(1/L*gp^2 + mu*(xs-xp)^2 - 2*mu/L*gp*(xp-xs));

% define polynomials
[s1, v1] = polynomial([fs fk fp xs xk gk gp],0);
[s2, v2] = polynomial([fs fk fp xs xk gk gp],0);
[s3, v3] = polynomial([fs fk fp xs xk gk gp],0);
[s4, v4] = polynomial([fs fk fp xs xk gk gp],0);
[s5, v5] = polynomial([fs fk fp xs xk gk gp],0);
[s6, v6] = polynomial([fs fk fp xs xk gk gp],0);

% feasible set
F = [sos(obj-[s1 s2 s3 s4 s5 s6]*[h1; h2; h3; h4; h5; h6]), ...
    sos(s1), sos(s2),sos(s3),sos(s4),sos(s5),sos(s6)];
[sol,v,Q] = solvesos(F,t,[],[v1;v2;v3;v4;v5;v6;t]);