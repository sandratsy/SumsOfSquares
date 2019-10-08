% Gradient Descent with Exact Line Search
% Searching for a degree-1 certificate using CVX
mu = 3; L = 5;
alpha = L/(2*(L-mu));

cvx_begin
    variable Q(9,9) semidefinite
    variable sig1 nonnegative
    variable sig2 nonnegative
    variable sig3 nonnegative
    variable sig4 nonnegative
    variable sig5 nonnegative
    variable sig6 nonnegative
    variables t the1 the2
    minimize(t)
    subject to
        t <= 1
        t >= 0
        0 == Q(1,1)
        1-t == 2*Q(1,2) - sig2 - sig4 + sig5 + sig6
        t == 2*Q(1,3) + sig1 + sig2 - sig3 - sig5
        -1 == 2*Q(1,4) - sig1 + sig3 + sig4 - sig6
        0 == Q(1,5:9)
        0 == Q(2:4,2:9)
        0 == Q(5,5) - alpha*mu*(sig2 + sig4 + sig5 + sig6)
        0 == 2*Q(5,6) + 2*alpha*mu*(sig2 + sig5)
        0 == 2*Q(5,7) + 2*alpha*mu*(sig4 + sig6)
        0 == 2*Q(5,8) - 2*alpha*mu/L*(sig2 + sig5) - sig5
        0 == 2*Q(5,9) - 2*alpha*mu/L*(sig4 + sig6) - sig6
        0 == Q(6,6) - alpha*mu*(sig1 + sig2 + sig3 + sig5)
        0 == 2*Q(6,7) + 2*alpha*mu*(sig1 + sig3)
        0 == 2*Q(6,8) + 2*alpha*mu/L*(sig1 + sig2 + sig3 + sig5) + sig3 + sig5
        0 == 2*Q(6,9) - sig1 - the1 - 2*alpha*mu/L*(sig1 + sig3)
        0 == Q(7,7) - alpha*mu*(sig1 + sig3 + sig4 + sig6)
        0 == 2*Q(7,8) - 2*alpha*mu/L*(sig1 + sig3) - sig3
        0 == 2*Q(7,9) + 2*alpha*mu/L*(sig1 + sig3 + sig4 + sig6) + sig1 + sig6 + the1
        0 == Q(8,8) - alpha/L*(sig1 + sig2 + sig3 + sig5)
        0 == 2*Q(8,9) + 2*alpha/L*(sig1 + sig3) + the2
        0 == Q(9,9) - alpha/L*(sig1 + sig3 + sig4 + sig6)
cvx_end