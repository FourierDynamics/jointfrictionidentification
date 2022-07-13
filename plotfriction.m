function plotfriction(optpara,qdot,tau_a)
a = optpara(1);
b = optpara(2);
S = optpara(3);
alpha = optpara(4);
v = optpara(5);
[m, n] = size(qdot);
tau = [];
for i= 1:m
    t = a*qdot(i) + b +S/(1 + exp(-alpha*(qdot(i)+v)));
    tau = [tau;t];
end
plot(qdot,tau_a);
hold on
plot(qdot,tau,'*r');