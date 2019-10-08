% Proximal Gradient Method with exact line search
% Searching for a degree-1 certificate using YALMIP
mu = 3; L = 5;
alpha = L/(2*(L-mu));

sdpvar fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp t
obj = t*(fk + hk) + (1-t)*(fs + hs) - (fp + hp);

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

% equality constraints (v == 0)
v13 = gp*gk + gp*sp + sp*gk + sp^2;
v14 = gp*xp - gp*xk + sp*xp - sp*xk;
v15 = gs^2 + 2*gs*ss + ss^2;

% define polynomials
[s1, c1] = polynomial([fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp],0);
[s2, c2] = polynomial([fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp],0);
[s3, c3] = polynomial([fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp],0);
[s4, c4] = polynomial([fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp],0);
[s5, c5] = polynomial([fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp],0);
[s6, c6] = polynomial([fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp],0);
[s7, c7] = polynomial([fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp],0);
[s8, c8] = polynomial([fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp],0);
[s9, c9] = polynomial([fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp],0);
[s10, c10] = polynomial([fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp],0);
[s11, c11] = polynomial([fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp],0);
[s12, c12] = polynomial([fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp],0);
[s13, c13] = polynomial([fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp],0);
[s14, c14] = polynomial([fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp],0);
[s15, c15] = polynomial([fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp],0);

% feasible set
F = [sos(obj-[s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13 s14 s15]*...
    [h1; h2; h3; h4; h5; h6; h7; h8; h9; h10; h11; h12; v13; v14; v15]), ...
    sos(s1), sos(s2),sos(s3),sos(s4),sos(s5),sos(s6),sos(s7),sos(s8),...
    sos(s9),sos(s10),sos(s11),sos(s12)];
[sol,v,Q] = solvesos(F,t,[],[c1;c2;c3;c4;c5;c6;c7;c8;c9;c10;c11;c12;c13;c14;c15;t]);