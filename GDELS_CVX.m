% Gradient Descent with Exact Line Search
% Searching for a mixed 0/1-degree certificate using CVX
mu = 3; L = 5;
alpha = L/(2*(L-mu));

cvx_begin
    variable Q(9,9) semidefinite
    variable l1 nonnegative
    variable l2 nonnegative
    variable l3 nonnegative
    variable l4 nonnegative
    variable l5 nonnegative
    variable l6 nonnegative
    variables l7 l8 t
    minimize(t)
    subject to
        t <= 1
        t >= 0
        0 == Q(1,1)
        1-t == 2*Q(1,2) - l2 - l4 + l5 + l6
        t == 2*Q(1,3) + l1 + l2 - l3 - l5
        -1 == 2*Q(1,4) - l1 + l3 + l4 - l6
        0 == Q(1,5:9)
        0 == Q(2:4,2:9)
        0 == Q(5,5) - alpha*mu*(l2 + l4 + l5 + l6)
        0 == 2*Q(5,6) + 2*alpha*mu*(l2 + l5)
        0 == 2*Q(5,7) + 2*alpha*mu*(l4 + l6)
        0 == 2*Q(5,8) - 2*alpha*mu/L*(l2 + l5) - l5
        0 == 2*Q(5,9) - 2*alpha*mu/L*(l4 + l6) - l6
        0 == Q(6,6) - alpha*mu*(l1 + l2 + l3 + l5)
        0 == 2*Q(6,7) + 2*alpha*mu*(l1 + l3)
        0 == 2*Q(6,8) + 2*alpha*mu/L*(l1 + l2 + l3 + l5) + l3 + l5
        0 == 2*Q(6,9) - l1 - l7 - 2*alpha*mu/L*(l1 + l3)
        0 == Q(7,7) - alpha*mu*(l1 + l3 + l4 + l6)
        0 == 2*Q(7,8) - 2*alpha*mu/L*(l1 + l3) - l3
        0 == 2*Q(7,9) + 2*alpha*mu/L*(l1 + l3 + l4 + l6) + l1 + l6 + l7
        0 == Q(8,8) - alpha/L*(l1 + l2 + l3 + l5)
        0 == 2*Q(8,9) + 2*alpha/L*(l1 + l3) + l8
        0 == Q(9,9) - alpha/L*(l1 + l3 + l4 + l6)
cvx_end