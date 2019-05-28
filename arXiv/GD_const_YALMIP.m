% Gradient Descent with Constant Step Size
% Searching for a degree-1 certificate using YALMIP

t = 0.3025; % adjust this value to desired degree of accuracy
mu = 3; L = 5;
alpha = L/(2*(L-mu));
gamma = 0.15; % set to any value in (0, 2/L)

% Updated variables
sdpvar fs fk fp xs xk gk gp
xp = xk - gamma*gk; % we substitute x_{k+1} with this equation

% p(z)
obj = t*(xk - xs)^2 - (xp - xs)^2;

% inequality constraints (c >= 0)
h1 = fk - fp - gp*(xk-xp) - alpha*(1/L*(gk-gp)^2 + mu*(xk-xp)^2 - 2*mu/L*(gp-gk)*(xp-xk));
h2 = fk - fs - alpha*(1/L*gk^2 + mu*(xk-xs)^2 + 2*mu/L*gk*(xs-xk));
h3 = fp - fk - gk*(xp-xk) - alpha*(1/L*(gp-gk)^2 + mu*(xp-xk)^2 - 2*mu/L*(gk-gp)*(xk-xp));
h4 = fp - fs - alpha*(1/L*gp^2 + mu*(xp-xs)^2 + 2*mu/L*gp*(xs-xp));
h5 = fs - fk - gk*(xs-xk) - alpha*(1/L*gk^2 + mu*(xs-xk)^2 - 2*mu/L*gk*(xk-xs));
h6 = fs - fp - gp*(xs-xp) - alpha*(1/L*gp^2 + mu*(xs-xp)^2 - 2*mu/L*gp*(xp-xs));

% feasible set
F = [h1>=0,h2>=0,h3>=0,h4>=0,h5>=0,h6>=0];

[info, X, moment, sos] = solvemoment(F,obj,[],1);

% if sos.t is numerically close to 0, the current t is feasible
sos.t