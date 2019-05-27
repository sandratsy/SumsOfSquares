% Proximal Gradient Method with constant step size
% Searching for a fully degree-1 certificate using YALMIP

t = 0.0625; % adjust this value to desired degree of accuracy
mu = 3; L = 5;
gamma = 2/(L+mu);
alpha = L/(2*(L-mu));

sdpvar fs fk fp hs hk hp xs xk gs gk gp sk sp
obj = t*(fk + hk) + (1-t)*(fs + hs) - (fp + hp);

% Eliminating variables
xp = xk - gamma*(gk + sp);
ss = -gs;

% inequality constraints (c >= 0)
h1 = fk - fp - gp*(xk-xp) - alpha*(1/L*(gk-gp)^2 + mu*(xk-xp)^2 - 2*mu/L*(gp-gk)*(xp-xk));
h2 = fk - fs - gs*(xk-xs) - alpha*(1/L*(gk-gs)^2 + mu*(xk-xs)^2 - 2*mu/L*(gs-gk)*(xs-xk));
h3 = fp - fk - gk*(xp-xk) - alpha*(1/L*(gp-gk)^2 + mu*(xp-xk)^2 - 2*mu/L*(gk-gp)*(xk-xp));
h4 = fp - fs - gs*(xp-xs) - alpha*(1/L*(gp-gs)^2 + mu*(xp-xs)^2 - 2*mu/L*(gs-gp)*(xs-xp));
h5 = fs - fk - gk*(xs-xk) - alpha*(1/L*(gs-gk)^2 + mu*(xs-xk)^2 - 2*mu/L*(gk-gs)*(xk-xs));
h6 = fs - fp - gp*(xs-xp) - alpha*(1/L*(gs-gp)^2 + mu*(xs-xp)^2 - 2*mu/L*(gp-gs)*(xp-xs));
h7 = hk - hp - sp*(xk - xp);
h8 = hk - hs - ss*(xk - xs);
h9 = hp - hk - sk*(xp - xk);
h10 = hp - hs - ss*(xp - xs);
h11 = hs - hk - sk*(xs - xk);
h12 = hs - hp - sp*(xs - xp);

% feasible set
F = [h1>=0,h2>=0,h3>=0,h4>=0,h5>=0,h6>=0,h7>=0,h8>=0,h9>=0,...
    h10>=0,h11>=0,h12>=0];

[info, X, moment, sos] = solvemoment(F,obj,[],1);

% if sos.t is numerically close to 0, the current t is feasible
sos.t

% solutions are found in sos.Q0, sos.Qi, etc.
