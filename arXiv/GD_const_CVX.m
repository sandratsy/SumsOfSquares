% Gradient Descent with Constant Step Size
% Searching for a degree-1 certificate using CVX
mu = 3; L = 5;
alpha = L/(2*(L-mu));
gamma = 0.4; % set to any value in (0, 2/L)

cvx_begin
    variable Q(8,8) semidefinite
    variable sig1 nonnegative
    variable sig2 nonnegative
    variable sig3 nonnegative
    variable sig4 nonnegative
    variable sig5 nonnegative
    variable sig6 nonnegative
    variables t
    minimize(t)
    subject to
        t <= 1
        t >= 0
        0 == Q(1,1)
        0 == 2*Q(1,2) - sig2 - sig4 + sig5 + sig6
        0 == 2*Q(1,3) + sig1 + sig2 - sig3 - sig5
        0 == 2*Q(1,4) - sig1 + sig3 + sig4 - sig6
        0 == Q(1,5:8)
        0 == Q(2:4,2:8)
        t-1 == Q(5,5) - alpha*mu*(sig2 + sig4 + sig5 + sig6)
        -2*t+2 == 2*Q(5,6) + 2*alpha*mu*(sig2 + sig4 + sig5 + sig6)
        -2*gamma == 2*Q(5,7) - 2*alpha*mu*((1/L)*(sig2 + sig5) + gamma*(sig4 + sig6)) - sig5
        0 == 2*Q(5,8) - 2*alpha*mu/L*(sig4 + sig6) - sig6
        t-1 == Q(6,6) - alpha*mu*(sig2 + sig4 + sig5 + sig6)
        2*gamma == 2*Q(6,7) + 2*alpha*mu*((1/L)*(sig2 + sig5) + gamma*(sig4 + sig6)) + sig5
        0 == 2*Q(6,8) + 2*alpha*mu/L*(sig4 + sig6) + sig6
        -gamma^2 == Q(7,7) - alpha/L*(sig1 + sig2 + sig3 + sig5) + gamma*sig3 ...
            - alpha*mu*gamma^2*(sig1 + sig3 + sig4 + sig6) + 2*alpha*gamma*mu/L*(sig1 + sig3)
        0 == 2*Q(7,8) - gamma*(sig1 + sig6) + 2*alpha/L*(sig1 + sig3 - gamma*mu*(sig1 + sig3 + sig4 + sig6))
        0 == Q(8,8) - alpha/L*(sig1 + sig3 + sig4 + sig6)
cvx_end

% Second SDP used to find sparse solutions 
cvx_begin
    variable Q(8,8) semidefinite
    variable sig1 nonnegative
    variable sig2 nonnegative
    variable sig3 nonnegative
    variable sig4 nonnegative
    variable sig5 nonnegative
    variable sig6 nonnegative
    minimize(sig1 + sig2 + sig3 + sig4 + sig5 + sig6)
    subject to
        0 == Q(1,1)
        0 == 2*Q(1,2) - sig2 - sig4 + sig5 + sig6
        0 == 2*Q(1,3) + sig1 + sig2 - sig3 - sig5
        0 == 2*Q(1,4) - sig1 + sig3 + sig4 - sig6
        0 == Q(1,5:8)
        0 == Q(2:4,2:8)
        t-1 == Q(5,5) - alpha*mu*(sig2 + sig4 + sig5 + sig6)
        -2*t+2 == 2*Q(5,6) + 2*alpha*mu*(sig2 + sig4 + sig5 + sig6)
        -2*gamma == 2*Q(5,7) - 2*alpha*mu*((1/L)*(sig2 + sig5) + gamma*(sig4 + sig6)) - sig5
        0 == 2*Q(5,8) - 2*alpha*mu/L*(sig4 + sig6) - sig6
        t-1 == Q(6,6) - alpha*mu*(sig2 + sig4 + sig5 + sig6)
        2*gamma == 2*Q(6,7) + 2*alpha*mu*((1/L)*(sig2 + sig5) + gamma*(sig4 + sig6)) + sig5
        0 == 2*Q(6,8) + 2*alpha*mu/L*(sig4 + sig6) + sig6
        -gamma^2 == Q(7,7) - alpha/L*(sig1 + sig2 + sig3 + sig5) + gamma*sig3 ...
            - alpha*mu*gamma^2*(sig1 + sig3 + sig4 + sig6) + 2*alpha*gamma*mu/L*(sig1 + sig3)
        0 == 2*Q(7,8) - gamma*(sig1 + sig6) + 2*alpha/L*(sig1 + sig3 - gamma*mu*(sig1 + sig3 + sig4 + sig6))
        0 == Q(8,8) - alpha/L*(sig1 + sig3 + sig4 + sig6)
cvx_end