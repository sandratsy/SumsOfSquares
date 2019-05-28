% Gradient Descent with Exact Line Search
% Verifying that the set of linear equalities hold for the optimal Q,
% sigmas, thetas and t

syms L mu
alpha = L/(2*(L-mu));

% optimal Q, t, sigmas and thetas
Q = zeros(9,9);
Q = sym(Q);
Q(5,5) = 2*L^2*mu^2/((L+mu)^2*(L-mu));
Q(5,6) = -L*mu^2/(L+mu)^2;
Q(5,7) = -L*mu^2/((L+mu)*(L-mu));
Q(5,8) = L*mu/(L+mu)^2;
Q(5,9) = L*mu/((L+mu)*(L-mu));
Q(6,6) = L*mu*(L + 3*mu)/(2*(L + mu)^2);
Q(6,7) = -L*mu/(2*(L + mu));
Q(6,8) = -mu*(3*L + mu)/(2*(L + mu)^2);
Q(6,9) = -mu/(2*(L + mu));
Q(7,7) = L*mu/(2*(L - mu));
Q(7,8) = mu/(2*(L + mu));
Q(7,9) = -mu/(2*(L - mu));
Q(8,8) = (L + 3*mu)/(2*(L + mu)^2);
Q(8,9) = 1/(2*(L + mu));
Q(9,9) = 1/(2*(L - mu));
for i = 5:9
    for j = i+1:9
        Q(j,i) = Q(i,j);
    end
end

t = ((L-mu)/(L+mu))^2;
sig1 = (L-mu)/(L+mu);
sig5 = 2*mu*(L-mu)/(L+mu)^2;
sig6 = 2*mu/(L+mu);
the1 = -1;
the2 = -2/(L+mu);

syms fs fk fp xs xk xp gk gp
zbar = [1; fs; fk; fp; xs; xk; xp; gk; gp];
zbarT = [1 fs fk fp xs xk xp gk gp];

% constraints (h >= 0, v == 0)
h1 = fk - fp - gp*(xk - xp) - alpha*(1/L*(gk-gp)^2 + mu*(xk-xp)^2 - 2*mu/L*(gp-gk)*(xp-xk));
h5 = fs - fk - gk*(xs - xk) - alpha*(1/L*gk^2 + mu*(xs-xk)^2 - 2*mu/L*gk*(xk-xs));
h6 = fs - fp - gp*(xs - xp) - alpha*(1/L*gp^2 + mu*(xs-xp)^2 - 2*mu/L*gp*(xp-xs));
v1 = gp*(xp - xk);
v2 = gp*gk;

% check whether linear equalities hold
pz = t*(fk - fs) - (fp - fs);
check = simplify(zbarT*Q*zbar + sig1*h1 + sig5*h5 + sig6*h6 + the1*v1 + the2*v2);
diff = simplify(pz - check); % should be zero