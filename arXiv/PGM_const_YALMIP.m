% Proximal Gradient Method with constant step size
% Searching for a degree-1 certificate using YALMIP
mu = 3; L = 5;
alpha = L/(2*(L-mu));
gamma = 2/(L+mu); % set to any value in (0, 2/L)

sdpvar fs fk fp hs hk hp xs xk gs gk gp sk sp t
obj = t*(gk + sk)^2 - (gp + sp)^2;

% Eliminating variables
xp = xk - gamma*(gk + sp);
ss = -gs;

% inequality constraints (h >= 0)
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

% define polynomials
[s1, v1] = polynomial([fs fk fp hs hk hp xs xk gs gk gp sk sp],0);
[s2, v2] = polynomial([fs fk fp hs hk hp xs xk gs gk gp sk sp],0);
[s3, v3] = polynomial([fs fk fp hs hk hp xs xk gs gk gp sk sp],0);
[s4, v4] = polynomial([fs fk fp hs hk hp xs xk gs gk gp sk sp],0);
[s5, v5] = polynomial([fs fk fp hs hk hp xs xk gs gk gp sk sp],0);
[s6, v6] = polynomial([fs fk fp hs hk hp xs xk gs gk gp sk sp],0);
[s7, v7] = polynomial([fs fk fp hs hk hp xs xk gs gk gp sk sp],0);
[s8, v8] = polynomial([fs fk fp hs hk hp xs xk gs gk gp sk sp],0);
[s9, v9] = polynomial([fs fk fp hs hk hp xs xk gs gk gp sk sp],0);
[s10, v10] = polynomial([fs fk fp hs hk hp xs xk gs gk gp sk sp],0);
[s11, v11] = polynomial([fs fk fp hs hk hp xs xk gs gk gp sk sp],0);
[s12, v12] = polynomial([fs fk fp hs hk hp xs xk gs gk gp sk sp],0);

% feasible set
F = [sos(obj-[s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12]*[h1; h2; h3; h4; h5; ...
    h6; h7; h8; h9; h10; h11; h12]), sos(s1),sos(s2),sos(s3),sos(s4),sos(s5),...
    sos(s6),sos(s7),sos(s8),sos(s9),sos(s10),sos(s11),sos(s12)];
[sol,v,Q] = solvesos(F,t,[],[v1;v2;v3;v4;v5;v6;v7;v8;v9;v10;v11;v12;t]);