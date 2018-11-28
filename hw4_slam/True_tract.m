clear;
load('data_set');%加载dataset数据
%--------------------------- PLOT RESULTS --------------------------------
figure; hold on;
plot(GPSLon(1:560),GPSLat(1:560));   % 机器人真实位置
% GPSLon 为机器人所在位置的纵坐标序列
legend('机器人真实位置');
xlabel('x [meters]'); ylabel('y [meters]');
axis([-10 20 -25 20]);
%------------------------ end of PLOT RESULTS -----------------------------
