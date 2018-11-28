clear ; 
load('data_set');
load('beac_juan3.mat');   % beac_juan3.mat������Ϊ·��,����whos �Cfile�鿴���ļ��е�����
landmarks=estbeac;
clear estbeac; % true landmark positions (measured w/GPS sensor)
%--------------------------- PLOT RESULTS --------------------------------
figure; hold on;
plot(landmarks(:,1),landmarks(:,2),'*');    % ·����ʵλ��
legend('·����ʵλ��');
xlabel('x [meters]');
ylabel('y [meters]');
axis([-10 20 -25 20]);
%------------------------ end of PLOT RESULTS -----------------------------
