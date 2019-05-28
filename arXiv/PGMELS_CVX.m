% Proximal Gradient Method with exact line search
% Searching for a degree-1 certificate using CVX

mu = 3; L = 5;
alpha = L/(2*(L-mu));

cvx_begin
    variable Q(16,16) semidefinite
    variables t the1 the2 the3
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
        1-t == 2*Q(1,2) - sig2 - sig4 + sig5 + sig6
        t == 2*Q(1,3) + sig1 + sig2 - sig3 - sig5
        -1 == 2*Q(1,4) - sig1 + sig3 + sig4 - sig6
        1-t == 2*Q(1,5) - sig8 - sig10 + sig11 + sig12
        t == 2*Q(1,6) + sig7 + sig8 - sig9 - sig11
        -1 == 2*Q(1,7) - sig7 + sig9 + sig10 - sig12
        0 == Q(1,8:16)
        0 == Q(2:7,2:16)
        0 == Q(8,8) - alpha*mu*(sig2 + sig4 + sig5 + sig6)
        0 == 2*Q(8,9) + 2*alpha*mu*(sig2 + sig5)
        0 == 2*Q(8,10) + 2*alpha*mu*(sig4 + sig6)
        0 == 2*Q(8,11) + sig2 + sig4 + 2*alpha*mu/L*(sig2 + sig4 + sig5 + sig6)
        0 == 2*Q(8,12) - 2*alpha*mu/L*(sig2 + sig5) - sig5
        0 == 2*Q(8,13) - 2*alpha*mu/L*(sig4 + sig6) - sig6
        0 == 2*Q(8,14) + sig8 + sig10
        0 == 2*Q(8,15) - sig11
        0 == 2*Q(8,16) - sig12
        0 == Q(9,9) - alpha*mu*(sig1 + sig2 + sig3 + sig5)
        0 == 2*Q(9,10) + 2*alpha*mu*(sig1 + sig3)
        0 == 2*Q(9,11) - sig2 - 2*alpha*mu/L*(sig2 + sig5)
        0 == 2*Q(9,12) + 2*alpha*mu/L*(sig1 + sig2 + sig3 + sig5) + sig3 + sig5
        0 == 2*Q(9,13) - the2 - sig1 - 2*alpha*mu/L*(sig1 + sig3)
        0 == 2*Q(9,14) - sig8
        0 == 2*Q(9,15) + sig9 + sig11
        0 == 2*Q(9,16) - the2 - sig7
        0 == Q(10,10) - alpha*mu*(sig1 + sig3 + sig4 + sig6)
        0 == 2*Q(10,11) - sig4 - 2*alpha*mu/L*(sig4 + sig6)
        0 == 2*Q(10,12) - 2*alpha*mu/L*(sig1 + sig3) - sig3
        0 == 2*Q(10,13) + the2 + sig1 + 2*alpha*mu/L*(sig1 + sig3 + sig4 + sig6) + sig6
        0 == 2*Q(10,14) - sig10
        0 == 2*Q(10,15) - sig9
        0 == 2*Q(10,16) + the2 + sig7 + sig12
        0 == Q(11,11) - alpha/L*(sig2 + sig4 + sig5 + sig6) + the3
        0 == 2*Q(11,12) + 2*alpha/L*(sig2 + sig5)
        0 == 2*Q(11,13) + 2*alpha/L*(sig4 + sig6)
        0 == 2*Q(11,14) + 2*the3
        0 == Q(11,15:16)
        0 == Q(12,12) - alpha/L*(sig1 + sig2 + sig3 + sig5)
        0 == 2*Q(12,13) + the1 + 2*alpha/L*(sig1 + sig3)
        0 == Q(12,14:15)
        0 == 2*Q(12,16) + the1
        0 == Q(13,13) - alpha/L*(sig1 + sig3 + sig4 + sig6)
        0 == Q(13,14:15)
        0 == 2*Q(13,16) + the1
        0 == Q(14,14) + the3
        0 == Q(14:15,15:16)
        0 == Q(16,16) + the1
cvx_end