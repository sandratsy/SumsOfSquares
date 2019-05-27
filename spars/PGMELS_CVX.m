% Proximal Gradient Method with exact line search
% Searching for a mixed 0/1-degree certificate using CVX

mu = 3; L = 5;
alpha = L/(2*(L-mu));

cvx_begin
    variable Q(16,16) semidefinite
    variables t l1 l2 l3
    variable l4 nonnegative
    variable l5 nonnegative
    variable l6 nonnegative
    variable l7 nonnegative
    variable l8 nonnegative
    variable l9 nonnegative
    variable l10 nonnegative
    variable l11 nonnegative
    variable l12 nonnegative
    variable l13 nonnegative
    variable l14 nonnegative
    variable l15 nonnegative
    minimize(t)
    subject to
        t <= 1
        t >= 0
        0 == Q(1,1)
        1-t == 2*Q(1,2) - l5 - l7 + l8 + l9
        t == 2*Q(1,3) + l4 + l5 - l6 - l8
        -1 == 2*Q(1,4) - l4 + l6 + l7 - l9
        1-t == 2*Q(1,5) - l11 - l13 + l14 + l15
        t == 2*Q(1,6) + l10 + l11 - l12 - l14
        -1 == 2*Q(1,7) - l10 + l12 + l13 - l15
        0 == Q(1,8:10)
        0 == 2*Q(1,11) + l3
        0 == Q(1,12:13)
        0 == 2*Q(1,14) + l3
        0 == Q(1,15:16)
        0 == Q(2:7,2:16)
        0 == Q(8,8) - alpha*mu*(l5 + l7 + l8 + l9)
        0 == 2*Q(8,9) + 2*alpha*mu*(l5 + l8)
        0 == 2*Q(8,10) + 2*alpha*mu*(l7 + l9)
        0 == 2*Q(8,11) + l5 + l7 + 2*alpha*mu/L*(l5 + l7 + l8 + l9)
        0 == 2*Q(8,12) - 2*alpha*mu/L*(l5 + l8) - l8
        0 == 2*Q(8,13) - 2*alpha*mu/L*(l7 + l9) - l9
        0 == 2*Q(8,14) + l11 + l13
        0 == 2*Q(8,15) - l14
        0 == 2*Q(8,16) - l15
        0 == Q(9,9) - alpha*mu*(l4 + l5 + l6 + l8)
        0 == 2*Q(9,10) + 2*alpha*mu*(l4 + l6)
        0 == 2*Q(9,11) - l5 - 2*alpha*mu/L*(l5 + l8)
        0 == 2*Q(9,12) + 2*alpha*mu/L*(l4 + l5 + l6 + l8) + l6 + l8
        0 == 2*Q(9,13) - l2 - l4 - 2*alpha*mu/L*(l4 + l6)
        0 == 2*Q(9,14) - l11
        0 == 2*Q(9,15) + l12 + l14
        0 == 2*Q(9,16) - l2 - l10
        0 == Q(10,10) - alpha*mu*(l4 + l6 + l7 + l9)
        0 == 2*Q(10,11) - l7 - 2*alpha*mu/L*(l7 + l9)
        0 == 2*Q(10,12) - 2*alpha*mu/L*(l4 + l6) - l6
        0 == 2*Q(10,13) + l2 + l4 + 2*alpha*mu/L*(l4 + l6 + l7 + l9) + l9
        0 == 2*Q(10,14) - l13
        0 == 2*Q(10,15) - l12
        0 == 2*Q(10,16) + l2 + l10 + l15
        0 == Q(11,11) - alpha/L*(l5 + l7 + l8 + l9)
        0 == 2*Q(11,12) + 2*alpha/L*(l5 + l8)
        0 == 2*Q(11,13) + 2*alpha/L*(l7 + l9)
        0 == Q(11,14:16)
        0 == Q(12,12) - alpha/L*(l4 + l5 + l6 + l8)
        0 == 2*Q(12,13) + l1 + 2*alpha/L*(l4 + l6)
        0 == Q(12,14:15)
        0 == 2*Q(12,16) + l1
        0 == Q(13,13) - alpha/L*(l4 + l6 + l7 + l9)
        0 == Q(13,14:15)
        0 == 2*Q(13,16) + l1
        0 == Q(14:15,14:16)
        0 == Q(16,16) + l1
cvx_end
