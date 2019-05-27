% Gradient Descent with Exact Line Search
% Searching for a fully degree-1 certificate using YALMIP

t = 0.0625; % adjust this value to desired degree of accuracy
mu = 3; L = 5;
alpha = L/(2*(L-mu));

sdpvar fs fk fp xs xk xp gk gp
obj = t*(fk - fs) - (fp - fs);

% inequality constraints (h >= 0)
h1 = fk - fp - gp*(xk-xp) - alpha*(1/L*(gk-gp)^2 + mu*(xk-xp)^2 - 2*mu/L*(gp-gk)*(xp-xk));
h2 = fk - fs - alpha*(1/L*gk^2 + mu*(xk-xs)^2 + 2*mu/L*gk*(xs-xk));
h3 = fp - fk - gk*(xp-xk) - alpha*(1/L*(gp-gk)^2 + mu*(xp-xk)^2 - 2*mu/L*(gk-gp)*(xk-xp));
h4 = fp - fs - alpha*(1/L*gp^2 + mu*(xp-xs)^2 + 2*mu/L*gp*(xs-xp));
h5 = fs - fk - gk*(xs-xk) - alpha*(1/L*gk^2 + mu*(xs-xk)^2 - 2*mu/L*gk*(xk-xs));
h6 = fs - fp - gp*(xs-xp) - alpha*(1/L*gp^2 + mu*(xs-xp)^2 - 2*mu/L*gp*(xp-xs));

% equality constraints (g == 0)
g1 = gp*(xp - xk);
g2 = gp*gk;

% feasible set
F = [h1>=0,h2>=0,h3>=0,h4>=0,h5>=0,h6>=0,g1==0,g2==0];

[info, X, moment, sos] = solvemoment(F,obj,[],1);

% if sos.t is numerically close to 0, the current t is feasible
sos.t

% solutions are found in sos.Q0, sos.Qi, etc.
