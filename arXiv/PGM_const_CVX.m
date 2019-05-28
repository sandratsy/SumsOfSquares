% Proximal Gradient Method with constant step size
% Searching for a degre-1 certificate using CVX

mu = 3; L = 5;
gamma = 0.05; % set to any value in (0, 2/L)
alpha = L/(2*(L-mu));

cvx_begin
    variable Q(14,14) semidefinite
    variables t
    variable sig1 nonnegative
    variable sig2 nonnegative
    variable sig3 nonnegative
    variable sig4 nonnegative
    variable sig5 nonnegative
    variable sig6 nonnegative
    variable sig7 nonnegative
    variable sig8 nonnegative
    variable sig9 nonnegative
    variable sig10 nonnegative
    variable sig11 nonnegative
    variable sig12 nonnegative
    minimize(t)
    subject to
        t <= 1
        t >= 0
        0 == Q(1,1)
        0 == 2*Q(1,2) - sig2 - sig4 + sig5 + sig6
        0 == 2*Q(1,3) + sig1 + sig2 - sig3 - sig5
        0 == 2*Q(1,4) - sig1 + sig3 + sig4 - sig6
        0 == 2*Q(1,5) - sig8 - sig10 + sig11 + sig12
        0 == 2*Q(1,6) + sig7 + sig8 - sig9 - sig11
        0 == 2*Q(1,7) - sig7 + sig9 + sig10 - sig12
        0 == Q(1,8:14)
        0 == Q(2:7,2:14)
        0 == Q(8,8) - alpha*mu*(sig2 + sig4 + sig5 + sig6)
        0 == 2*Q(8,9) + 2*alpha*mu*(sig2 + sig4 + sig5 + sig6)
        0 == 2*Q(8,10) + sig2 + sig4 - sig8 - sig10 + 2*alpha*mu/L*(sig2 + sig4 + sig5 + sig6)
        0 == 2*Q(8,11) - 2*alpha*mu*((1/L)*(sig2 + sig5) + gamma*(sig4 + sig6)) - sig5
        0 == 2*Q(8,12) - 2*alpha*mu/L*(sig4 + sig6) - sig6
        0 == 2*Q(8,13) - sig11
        0 == 2*Q(8,14) - 2*alpha*mu*gamma*(sig4 + sig6) - sig12
        0 == Q(9,9) - alpha*mu*(sig2 + sig4 + sig5 + sig6)
        0 == 2*Q(9,10) - sig2 - sig4 - 2*alpha*mu/L*(sig2 + sig4 + sig5 + sig6) + sig8 + sig10
        0 == 2*Q(9,11) + 2*alpha*mu*((1/L)*(sig2 + sig5) + gamma*(sig4 + sig6)) + sig5
        0 == 2*Q(9,12) + 2*alpha*mu/L*(sig4 + sig6) + sig6
        0 == 2*Q(9,13) + sig11
        0 == 2*Q(9,14) + 2*alpha*mu*gamma*(sig4 + sig6) + sig12
        0 == Q(10,10) - alpha/L*(sig2 + sig4 + sig5 + sig6)
        0 == 2*Q(10,11) + 2*alpha/L*(sig2 + sig5 + mu*gamma*(sig4 + sig6)) + gamma*(sig4 - sig10)
        0 == 2*Q(10,12) + 2*alpha/L*(sig4 + sig6)
        0 == Q(10,13)
        0 == 2*Q(10,14) + gamma*(sig4 - sig10) + 2*alpha*mu/L*gamma*(sig4 + sig6)
        t == Q(11,11) - alpha/L*(sig1 + sig2 + sig3 + sig5) - alpha*mu*gamma^2*(sig1 + sig3 + sig4 + sig6) ...
            + 2*alpha*gamma*mu/L*(sig1 + sig3) + gamma*sig3
        0 == 2*Q(11,12) - gamma*(sig1 + sig6) + 2*alpha/L*(sig1 + sig3 - mu*gamma*(sig1 + sig3 + sig4 + sig6))
        2*t == 2*Q(11,13) + gamma*sig9
        0 == 2*Q(11,14) + 2*alpha*mu/L*gamma*(sig1 + sig3) + gamma*(sig3 - sig7 - sig12) ...
            - 2*alpha*mu*gamma^2*(sig1 + sig3 + sig4 + sig6)
        -1 == Q(12,12) - alpha/L*(sig1 + sig3 + sig4 + sig6)
        0 == Q(12,13)
        -2 == 2*Q(12,14) - gamma*(sig1 + sig6) - 2*alpha*mu/L*gamma*(sig1 + sig3 + sig4 + sig6)
        t == Q(13,13)
        0 == 2*Q(13,14) + gamma*sig9
        -1 == Q(14,14) - alpha*mu*gamma^2*(sig1 + sig3 + sig4 + sig6) - gamma*(sig7 + sig12)
cvx_end