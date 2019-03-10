% Proximal Gradient Method with constant step size
% Searching for a mixed 0/1-degree certificate using CVX

mu = 3; L = 5;
gamma = 2/(L+mu);
alpha = L/(2*(L-mu));

cvx_begin
    variable Q(14,14) semidefinite
    variables t
    variable l1 nonnegative
    variable l2 nonnegative
    variable l3 nonnegative
    variable l4 nonnegative
    variable l5 nonnegative
    variable l6 nonnegative
    variable l7 nonnegative
    variable l8 nonnegative
    variable l9 nonnegative
    variable l10 nonnegative
    variable l11 nonnegative
    variable l12 nonnegative
    minimize(t)
    subject to
        t <= 1
        t >= 0
        0 == Q(1,1)
        1-t == 2*Q(1,2) - l2 - l4 + l5 + l6
        t == 2*Q(1,3) + l1 + l2 - l3 - l5
        -1 == 2*Q(1,4) - l1 + l3 + l4 - l6
        1-t == 2*Q(1,5) - l8 - l10 + l11 + l12
        t == 2*Q(1,6) + l7 + l8 - l9 - l11
        -1 == 2*Q(1,7) - l7 + l9 + l10 - l12
        0 == Q(1,8:14)
        0 == Q(2:7,2:14)
        0 == Q(8,8) - alpha*mu*(l2 + l4 + l5 + l6)
        0 == 2*Q(8,9) + 2*alpha*mu*(l2 + l4 + l5 + l6)
        0 == 2*Q(8,10) + l2 + l4 - l8 - l10 + 2*alpha*mu/L*(l2 + l4 + l5 + l6)
        0 == 2*Q(8,11) - 2*alpha*mu*((1/L)*(l2 + l5) + gamma*(l4 + l6)) - l5
        0 == 2*Q(8,12) - 2*alpha*mu/L*(l4 + l6) - l6
        0 == 2*Q(8,13) - l11
        0 == 2*Q(8,14) - 2*alpha*mu*gamma*(l4 + l6) - l12
        0 == Q(9,9) - alpha*mu*(l2 + l4 + l5 + l6)
        0 == 2*Q(9,10) - l2 - l4 - 2*alpha*mu/L*(l2 + l4 + l5 + l6) + l8 + l10
        0 == 2*Q(9,11) + 2*alpha*mu*((1/L)*(l2 + l5) + gamma*(l4 + l6)) + l5
        0 == 2*Q(9,12) + 2*alpha*mu/L*(l4 + l6) + l6
        0 == 2*Q(9,13) + l11
        0 == 2*Q(9,14) + 2*alpha*mu*gamma*(l4 + l6) + l12
        0 == Q(10,10) - alpha/L*(l2 + l4 + l5 + l6)
        0 == 2*Q(10,11) + 2*alpha/L*(l2 + l5 + mu*gamma*(l4 + l6)) + gamma*(l4 - l10)
        0 == 2*Q(10,12) + 2*alpha/L*(l4 + l6)
        0 == Q(10,13)
        0 == 2*Q(10,14) + gamma*(l4 - l10) + 2*alpha*mu/L*gamma*(l4 + l6)
        0 == Q(11,11) - alpha/L*(l1 + l2 + l3 + l5) - alpha*mu*gamma^2*(l1 + l3 + l4 + l6) ...
            + 2*alpha*gamma*mu/L*(l1 + l3) + gamma*l3
        0 == 2*Q(11,12) - gamma*(l1 + l6) + 2*alpha/L*(l1 + l3 - mu*gamma*(l1 + l3 + l4 + l6))
        0 == 2*Q(11,13) + gamma*l9
        0 == 2*Q(11,14) + 2*alpha*mu/L*gamma*(l1 + l3) + gamma*(l3 - l7 - l12) ...
            - 2*alpha*mu*gamma^2*(l1 + l3 + l4 + l6)
        0 == Q(12,12) - alpha/L*(l1 + l3 + l4 + l6)
        0 == Q(12,13)
        0 == 2*Q(12,14) - gamma*(l1 + l6) - 2*alpha*mu/L*gamma*(l1 + l3 + l4 + l6)
        0 == Q(13,13)
        0 == 2*Q(13,14) + gamma*l9
        0 == Q(14,14) - alpha*mu*gamma^2*(l1 + l3 + l4 + l6) - gamma*(l7 + l12)
cvx_end