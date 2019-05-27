% Proximal Gradient Method with constant step size
% Verifying that the set of linear equalities hold for the optimal Q,
% lambdas and t

syms L mu
alpha = L/(2*(L-mu));
gamma = 2/(L+mu);

% optimal Q, t and lambdas
Q = zeros(14,14);
Q = sym(Q);
Q(8,8) = 2*L^2*mu^2/((L + mu)^2*(L - mu));
Q(8,9) = -2*L^2*mu^2/((L + mu)^2*(L - mu));
Q(8,10) = -2*L*mu^2/((L + mu)^2*(L - mu));
Q(8,11) = L*mu/((L + mu)*(L - mu));
Q(8,12) = L*mu/((L + mu)*(L - mu));
Q(8,14) = 2*L^2*mu/((L + mu)^2*(L - mu));
Q(9,9) = 2*L^2*mu^2/((L + mu)^2*(L - mu));
Q(9,10) = 2*L*mu^2/((L + mu)^2*(L - mu));
Q(9,11) = -L*mu/((L+mu)*(L-mu));
Q(9,12) = -L*mu/((L+mu)*(L-mu));
Q(9,14) = -2*L^2*mu/((L + mu)^2*(L - mu));
Q(10,10) = 2*L*mu/((L + mu)^2*(L - mu));
Q(10,11) = -mu/((L+mu)*(L-mu));
Q(10,12) = -mu/((L+mu)*(L-mu));
Q(10,14) = -2*mu^2/((L + mu)^2*(L - mu));
Q(11,11) = 1/(2*(L - mu));
Q(11,12) = 1/(2*(L - mu));
Q(11,14) = L/((L+mu)*(L-mu));
Q(12,12) = 1/(2*(L - mu));
Q(12,14) = L/((L+mu)*(L-mu));
Q(14,14) = 2*(L^2 + L*mu - mu^2)/((L + mu)^2*(L - mu));
for i = 8:14
    for j = i+1:14
        Q(j,i) = Q(i,j);
    end
end

t = ((L-mu)/(L+mu))^2;
l1 = (L-mu)/(L+mu);
l5 = 2*mu*(L-mu)/(L+mu)^2;
l6 = 2*mu/(L+mu);
l7 = ((L-mu)/(L+mu))^2;
l12 = 4*L*mu/(L+mu)^2;

syms fs fk fp hs hk hp xs xk gs gk gp sk sp;
zbar = [1; fs; fk; fp; hs; hk; hp; xs; xk; gs; gk; gp; sk; sp];
zbarT = [1 fs fk fp hs hk hp xs xk gs gk gp sk sp];

% Eliminating variables using equality constraints
xp = xk - gamma*(gk + sp);
ss = -gs;

% inequality constraints (h >= 0)
h1 = fk - fp - gp*(xk-xp) - alpha*(1/L*(gk-gp)^2 + mu*(xk-xp)^2 - 2*mu/L*(gp-gk)*(xp-xk));
h5 = fs - fk - gk*(xs-xk) - alpha*(1/L*(gs-gk)^2 + mu*(xs-xk)^2 - 2*mu/L*(gk-gs)*(xk-xs));
h6 = fs - fp - gp*(xs-xp) - alpha*(1/L*(gs-gp)^2 + mu*(xs-xp)^2 - 2*mu/L*(gp-gs)*(xp-xs));
h7 = hk - hp - sp*(xk - xp);
h12 = hs - hp - sp*(xs - xp);

% check whether linear equalities hold
pz = t*(fk + hk) + (1-t)*(fs + hs) - (fp + hp);
check = simplify(zbarT*Q*zbar + l1*h1 + l5*h5 + l6*h6 + l7*h7 + l12*h12);
diff = simplify(pz - check); % should be zero
