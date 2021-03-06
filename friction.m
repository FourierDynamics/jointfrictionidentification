
a1 = load("datacollection_5_0.05_plus.txt");
a2 = load("datacollection_5_0.05_minus.txt");
ki = 0.067;
kr = 51;

fc = 5;
fs = 250;
filter5=designfilt('lowpassfir', 'PassbandFrequency', 0.8*fc, 'StopbandFrequency', fc, 'PassbandRipple', 0.01, 'StopbandAttenuation', 60, 'SampleRate', fs);
% motor15_qVelEst_1_001_plus_filterfir5=filtfilt(filter5,motor15_qVelEst_1_001_plus);
% ip = [v i]
ip = [a1(:,14) a1(:,17);a2(:,14) a2(:,17);];
figure;
plot(ip(:,1),ip(:,2)*ki*kr);
hold on;
ip_fft_v = filtfilt(filter5,ip(:,1));
ip_fft_i = filtfilt(filter5,ip(:,2));

plot(ip_fft_v,ip_fft_i*ki*kr,'-r');

%%
    options = optimoptions('fmincon', 'PlotFcn','optimplotfval');   %画目标函数值
    options.MaxFunctionEvaluations = 1000000;   %最大目标函数计算次数
    options.MaxIterations =10000;%最大迭代次数
    options.Display = 'iter';
    options.StepTolerance = 1e-8;
    options.OptimalityTolerance = 1e-7;
    options.FunctionTolerance = 1e-7;

init_para = rand(1)*ones(5,1);
init_para'

A= [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = [];
[Opt_para, Opt_obj, flag_exit] = fmincon(@(qp_para) fun_t(qp_para, ip_fft_v,ip_fft_i*ki*kr),...
                                                    init_para, A, b, Aeq, beq, lb, ub, ...
                                                    nonlcon, options); 
                                               
%%
% 
plotfriction(Opt_para,ip_fft_v,ip_fft_i*ki*kr);

                                                
                                                
                                                
                                                
                                                
                                                
                                                
                                                
                                                