% Gradient Descent with Armijo Rule
% Verifying that the set of linear equalities hold for the optimal Q,
% sigmas and t

syms L mu eps eta 
alpha = L/(2*(L-mu));

% optimal Q, t, sigmas
Q = zeros(9,9);
Q = sym(Q);
Q(5,5) = 2*mu^2*eps*(1-eps)/(eta*(L-mu));
Q(5,6) = -2*mu^2*eps*(1-eps)/(eta*(L-mu));
Q(5,8) = 2*mu*eps*(1-eps)/(eta*(L-mu));
Q(6,6) = 2*mu^2*eps*(1-eps)/(eta*(L-mu));
Q(6,8) = -2*mu*eps*(1-eps)/(eta*(L-mu));
Q(8,8) = 2*eps*(1-eps)/(eta*(L-mu));
for i = 5:9
    for j = i+1:9
        Q(j,i) = Q(i,j);
    end
end

t = 1 - 4*mu*eps*(1-eps)/(eta*L);
sig5 = 4*mu*eps*(1-eps)/(eta*L);
sig7 = 1;

syms fs fk fp xs xk xp gk gp
zbar = [1; fs; fk; fp; xs; xk; xp; gk; gp];
zbarT = [1 fs fk fp xs xk xp gk gp];

% constraints (h >= 0)
h5 = fs - fk - gk*(xs - xk) - alpha*(1/L*gk^2 + mu*(xs-xk)^2 - 2*mu/L*gk*(xk-xs));
h7 = fk - fp - 2*eps*(1-eps)/(eta*L)*gk^2;

% check whether linear equalities hold
pz = t*(fk - fs) - (fp - fs);
check = simplify(zbarT*Q*zbar + sig5*h5 + sig7*h7);
diff = simplify(pz - check); % should be zero