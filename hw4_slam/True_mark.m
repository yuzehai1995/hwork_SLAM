clear ; 
load('data_set');
load('beac_juan3.mat');   % beac_juan3.mat中内容为路标,命令whos –file查看该文件中的内容
landmarks=estbeac;
clear estbeac; % true landmark positions (measured w/GPS sensor)
%--------------------------- PLOT RESULTS --------------------------------
figure; hold on;
plot(landmarks(:,1),landmarks(:,2),'*');    % 路标真实位置
legend('路标真实位置');
xlabel('x [meters]');
ylabel('y [meters]');
axis([-10 20 -25 20]);
%------------------------ end of PLOT RESULTS -----------------------------
