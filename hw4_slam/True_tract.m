clear;
load('data_set');%����dataset����
%--------------------------- PLOT RESULTS --------------------------------
figure; hold on;
plot(GPSLon(1:560),GPSLat(1:560));   % ��������ʵλ��
% GPSLon Ϊ����������λ�õ�����������
legend('��������ʵλ��');
xlabel('x [meters]'); ylabel('y [meters]');
axis([-10 20 -25 20]);
%------------------------ end of PLOT RESULTS -----------------------------
