% Proximal Gradient Method with constant step size
% Verifying that equation (79) holds for the optimal Q, sigmas, thetas and t

syms L mu
alpha = L/(2*(L-mu));

% optimal Q, t and lambdas
Q = zeros(16,16);
Q = sym(Q);
Q(8,8) = 2*L^2*mu^2/((L+mu)^2 * (L-mu));
Q(8,9) = -L*mu^2/(L+mu)^2;
Q(8,10) = -L*mu^2/((L-mu)*(L+mu));
Q(8,11) = -2*L*mu^2/((L+mu)^2 * (L-mu));
Q(8,12) = L*mu/(L+mu)^2;
Q(8,13) = L*mu/((L-mu)*(L+mu));
Q(8,16) = 2*L*mu/(L+mu)^2;
Q(9,9) = L*mu*(L+3*mu)/(2*(L+mu)^2);
Q(9,10) = -L*mu/(2*(L+mu));
Q(9,11) = mu^2/(L+mu)^2;
Q(9,12) = -mu*(3*L+mu)/(2*(L+mu)^2);
Q(9,13) = -mu/(2*(L+mu));
Q(9,16) = -2*L*mu/(L+mu)^2;
Q(10,10) = L*mu/(2*(L-mu));
Q(10,11) = mu^2/((L-mu)*(L+mu));
Q(10,12) = mu/(2*(L+mu));
Q(10,13) = -mu/(2*(L-mu));
Q(11,11) = 2*L*mu/((L+mu)^2 * (L-mu));
Q(11,12) = -mu/(L+mu)^2;
Q(11,13) = -mu/((L-mu)*(L+mu));
Q(12,12) = (L+3*mu)/(2*(L+mu)^2);
Q(12,13) = 1/(2*(L+mu));
Q(12,16) = 1/(L+mu);
Q(13,13) = 1/(2*(L-mu));
Q(13,16) = 1/(L+mu);
Q(16,16) = 2/(L+mu);
for i = 8:16
    for j = i+1:16
        Q(j,i) = Q(i,j);
    end
end

t = ((L-mu)/(L+mu))^2;
the1 = -2/(L+mu);
the2 = -1;
sig1 = (L-mu)/(L+mu);
sig5 = 2*mu*(L-mu)/(L+mu)^2;
sig6 = 2*mu/(L+mu);
sig7 = ((L-mu)/(L+mu))^2;
sig12 = 4*L*mu/(L+mu)^2;

syms fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp
zbar = [1; fs; fk; fp; hs; hk; hp; xs; xk; xp; gs; gk; gp; ss; sk; sp];
zbarT = [1 fs fk fp hs hk hp xs xk xp gs gk gp ss sk sp];

% constraints (h >= 0, v == 0)
h1 = fk - fp - gp*(xk - xp) - alpha*(1/L*(gk-gp)^2 + mu*(xk-xp)^2 - 2*mu/L*(gp-gk)*(xp-xk));
h5 = fs - fk - gk*(xs - xk) - alpha*(1/L*(gs-gk)^2 + mu*(xs-xk)^2 - 2*mu/L*(gk-gs)*(xk-xs));
h6 = fs - fp - gp*(xs - xp) - alpha*(1/L*(gs-gp)^2 + mu*(xs-xp)^2 - 2*mu/L*(gp-gs)*(xp-xs));
h7 = hk - hp - sp*(xk - xp);
h12 = hs - hp - sp*(xs - xp);
v1 = gp*gk + gp*sp + sp*gk + sp^2;
v2 = gp*xp - gp*xk + sp*xp - sp*xk;

% check whether equation (79) holds
pz = t*(fk + hk) + (1-t)*(fs + hs) - (fp + hp);
check = simplify(zbarT*Q*zbar + the1*v1 + the2*v2 + sig1*h1 + sig5*h5 + sig6*h6 + sig7*h7 + sig12*h12);
diff = simplify(pz - check); %should be zero