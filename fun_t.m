function error = fun_t(optpara, q_dot_a, tau_a)
a = optpara(1);
b = optpara(2);
S = optpara(3);
alpha = optpara(4);
v = optpara(5);
[m , n] = size(q_dot_a);
tau = [];
for i=1:m
    t = a*q_dot_a(i) + b +S/(1 + exp(-alpha*(q_dot_a(i)+v)));
    tau = [tau;t];
end
error = norm(tau - tau_a);