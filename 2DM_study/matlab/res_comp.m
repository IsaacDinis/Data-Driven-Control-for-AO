res_cascaded_standalone = load('../data/res_cascaded_standalone.mat').res_cascaded_standalone;
res_cascaded_integrated = load('../data/res_cascaded_integrated.mat').res_cascaded_integrated;
res_1dm = load('../data/res_1dm.mat').res_1dm;

t = 0:1/4000:0.25-1/4000;
figure()
plot(t,res_cascaded_standalone(2,9001:10000))

hold on;
plot(t,res_cascaded_integrated(2,9001:10000),'LineWidth',1)



figure()
plot(t,res_cascaded_integrated(2,9001:10000),'LineWidth',1)
hold on;
plot(t,res_1dm(2,9001:10000))
% plot(t,res_1dm(2,9001:13000))