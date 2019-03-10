% Gradient Descent with Exact Line Search
% Verifying that the set of linear equalities hold for the optimal Q,
% lambdas and t

syms L mu
alpha = L/(2*(L-mu));

% optimal Q, t and lambdas
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
l1 = (L-mu)/(L+mu);
l5 = 2*mu*(L-mu)/(L+mu)^2;
l6 = 2*mu/(L+mu);
l7 = -1;
l8 = -2/(L+mu);

syms fs fk fp xs xk xp gk gp
zbar = [1; fs; fk; fp; xs; xk; xp; gk; gp];
zbarT = [1 fs fk fp xs xk xp gk gp];

% constraints
c1 = fk - fp - gp*(xk - xp) - alpha*(1/L*(gk-gp)^2 + mu*(xk-xp)^2 - 2*mu/L*(gp-gk)*(xp-xk));
c5 = fs - fk - gk*(xs - xk) - alpha*(1/L*gk^2 + mu*(xs-xk)^2 - 2*mu/L*gk*(xk-xs));
c6 = fs - fp - gp*(xs - xp) - alpha*(1/L*gp^2 + mu*(xs-xp)^2 - 2*mu/L*gp*(xp-xs));
c7 = gp*(xp - xk);
c8 = gp*gk;

% check whether linear equalities hold
pz = t*(fk - fs) - (fp - fs);
check = simplify(zbarT*Q*zbar + l1*c1 + l5*c5 + l6*c6 + l7*c7 + l8*c8);
diff = simplify(pz - check); % should be zero