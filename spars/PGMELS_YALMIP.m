% Proximal Gradient Method with exact line search
% Searching for a fully degree-1 certificate using YALMIP

t = 0.0625; % adjust this value to desired degree of accuracy
mu = 3; L = 5;
alpha = L/(2*(L-mu));

sdpvar fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp
obj = t*(fk + hk) + (1-t)*(fs + hs) - (fp + hp);

% equality constraints (g == 0)
g1 = gp*gk + gp*sp + sp*gk + sp^2;
g2 = gp*xp - gp*xk + sp*xp - sp*xk;
g3 = gs + ss;

% inequality constraints (h >= 0)
h4 = fk - fp - gp*(xk-xp) - alpha*(1/L*(gk-gp)^2 + mu*(xk-xp)^2 - 2*mu/L*(gp-gk)*(xp-xk));
h5 = fk - fs - gs*(xk-xs) - alpha*(1/L*(gk-gs)^2 + mu*(xk-xs)^2 - 2*mu/L*(gs-gk)*(xs-xk));
h6 = fp - fk - gk*(xp-xk) - alpha*(1/L*(gp-gk)^2 + mu*(xp-xk)^2 - 2*mu/L*(gk-gp)*(xk-xp));
h7 = fp - fs - gs*(xp-xs) - alpha*(1/L*(gp-gs)^2 + mu*(xp-xs)^2 - 2*mu/L*(gs-gp)*(xs-xp));
h8 = fs - fk - gk*(xs-xk) - alpha*(1/L*(gs-gk)^2 + mu*(xs-xk)^2 - 2*mu/L*(gk-gs)*(xk-xs));
h9 = fs - fp - gp*(xs-xp) - alpha*(1/L*(gs-gp)^2 + mu*(xs-xp)^2 - 2*mu/L*(gp-gs)*(xp-xs));
h10 = hk - hp - sp*(xk - xp);
h11 = hk - hs - ss*(xk - xs);
h12 = hp - hk - sk*(xp - xk);
h13 = hp - hs - ss*(xp - xs);
h14 = hs - hk - sk*(xs - xk);
h15 = hs - hp - sp*(xs - xp);

% feasible set
F = [g1==0, g2==0, g3==0,h4>=0,h5>=0,h6>=0,h7>=0,h8>=0,h9>=0,h10>=0,h11>=0,...
    h12>=0,h13>=0,h14>=0,h15>=0];

[info, X, moment, sos] = solvemoment(F,obj,[],1);

% if sos.t is numerically close to 0, the current t is feasible
sos.t

% solutions are found in sos.Q0, sos.Qi, etc.
