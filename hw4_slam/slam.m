%-------------------------------------------------------------------------
% FILE: slam.m 
% DESC: Implements the SLAM algorithm based on a robot model and sensor 
%       data from Nebot.
%landmarks(真实路标) 
%estbeac=estimation(估计)+beacon(路标)，将estbeac的值赋给landmarks然后将estbeac清除。
%xest                                   状态矩阵
%Pest                                   协方差矩阵
%A W Jh                                 雅可比矩阵
%Q                                      过程激励噪声协方差矩阵
%k                                      迭代次数
%K                                      卡尔曼增益
%numStates                              状态的数量  初始化时候就三个，即位置(x，y)和姿态(phi)
%phi                                    姿态，即与x轴的夹角
%Steering                               转向角
%alpha                                  转向角
%Ve = Velocity                          后轮的转速
%Vc                                     车水平平移速度
%zlaser=[rangeArray; bearingArray];     机器人与路标的距离和方位
%--------------------------------------------------------------------------
clear; clc;
load('data_set');
load('beac_juan3.mat');   %beac_juan3.mat中内容为路标,命令whos –file查看该文件中的内容
landmarks=estbeac; clear estbeac; % true landmark positions (measured w/GPS sensor)
%---------------------  DATA INITIALIZATION   -----------------------------
% Vehicle constants
L=2.83;     % distance between front and rear axles(前后轮轴之间的距离)
h=0.76;     % distance between center of rear axle and encoder(后轮轴中心与编码器的距离)
b=0.5;      % vertical distance from rear axle to laser(后轮轴中心到激光传感器的垂直距离)
a=3.78;     % horizontal distance from rear axle to laser(后轮轴中心到激光传感器的水平距离)

% EKF components used in multiple functions
global xest Pest A Q W   %定义全局变量
global numStates k  

xest = zeros(3,1);  % state matrix
Pest = zeros(3,3);  % covariance matrix2
numStates=3;        % initially, there are no landmarks(i.e. the number of value in the state matrix)
% matrixes initialization
A = zeros(3,3);                 % process model jacobian
W = eye(3);                     % process noise jacobian
Q = [.49   0   0; 0 .49 0; 0 0 (7*pi/180)^2];
Jh = zeros(2,3);                % sensor model jacobian
K = zeros(3,3);                 % kalman gain
 
% initialize time
To = min([TimeGPS(1),Time_VS(1),TimeLaser(1)]); % data starts being collected at 45 secs for some reason??
Time = Time - To;
t=Time(1);
tfinal = 112;   % there is 112 seconds worth of data
k=1;            % loop iteration counter(循环迭代计数)
 
% initialize estimates
xest(1,1) = GPSLon(1);      % x position of vehicle
xest(2,1) = GPSLat(1);      % y position of vehicle
xest(3,1) = -126*pi/180;    % phi of vehicle       % LB(?):how to get this value?
Pest(:,:) = [0.1, 0, 0; 0, 0.1, 0; 0, 0, 15*pi/180];
 
% initial inputs
alpha = Steering(1);    % measured steering angle(测量转向角)
Ve = Velocity(1);       % measured velocity at back wheel(测量后轮速率)
Vc = Velocity(1) / (1 - (h/L)*tan(alpha));  % measured velocity transformed to vehicle center
%--------------------- end of DATA INITIALIZATION  ------------------------
 
 
disp('Running through Kalman filter...');
 
 
n_dsp = 0;
%------------------  EXTENDED KALMAN FILTER PROCEDURE  --------------------
while t<tfinal
 
    k=k+1;
    t=Time(k);
    dt=t-Time(k-1);
    % Calculate time update matrices (prediction step of kalman filter)
    % Vc and alpha are updated after prediction...therefore, Vc is Vc(k-1) & alpha is alpha(k-1)
    
    %下面是公式（1）
    phi = xest(3,k-1);  % so I don't have to type in xest(3,k-1) for every phi below
    xest(1,k) = xest(1,k-1) + dt*Vc*cos(phi) - dt*Vc/L*tan(alpha)*(a*sin(phi)+b*cos(phi));     % x prediction
    xest(2,k) = xest(2,k-1) + dt*Vc*sin(phi) + dt*Vc/L*tan(alpha)*(a*cos(phi)-b*sin(phi));     % y prediction
    xest(3,k) = xest(3,k-1) + dt*Vc/L*tan(alpha);     % phi prediction
    xest(3,k) = normalizeAngle(xest(3,k));      % keep phi between -180 and 180 degrees
    
    % landmarks are assumed to be static, i.e. pi(k) = pi(k-1)         
    if(numStates>3)
        for i=4:numStates
            xest(i,k)= xest(i,k-1);
        end
    end
    
    % Calculate Jacobian A based on vehicle model (df/dx)
    A(1,1)=1; A(1,2)=0; A(1,3) = -dt*Vc*sin(phi) - dt*Vc/L*tan(alpha)*(a*cos(phi)-b*sin(phi));
    A(2,1)=0; A(2,2)=1; A(2,3) = dt*Vc*cos(phi) - dt*Vc/L*tan(alpha)*(a*sin(phi)+b*cos(phi)); 
    A(3,1)=0; A(3,2)=0; A(3,3)=1;
    % rest of A (which is just 2 null matrices and 1 identity) is added in new_state function
        
    % Calculate prediction for covariance matrix
    Pest = A*Pest*A' + W*Q*W';      %(公式2)
 
    % Check to see if new velocity and steering measurements are available
    if Sensor(k)==2   % encoders(编码器) are sensor #2
        Ve = Velocity(Index(k));
        alpha = Steering(Index(k));
        Vc = Ve / (1-(h/L)*tan(alpha));
    end  
    
    % Check to see if new laser(激光) data is available    
    if (Sensor(k)==3 && t>1)   % SICK is sensor #3                
        
        % from a single 180 degree scan, determine # of landmarks & center of each landmark
        [rangeArray ,bearingArray] = getLandmarkCenter(Laser(Index(k),:), Intensity(Index(k),:));   
        zlaser=[rangeArray; bearingArray];          % range/bearing data in observation model format
        numLandmarks = size(rangeArray,2);          % number of LMs captured in a single scan      
              
        for i=1:numLandmarks                            % for i = 1 to # of landmakrs picked up in single scan
            if(numStates==3)                            % if 1st observed landmark, update state and covariance
                new_state(zlaser(:,i));                 % add new state (i.e. increase all matrix sizes)
                numStates = numStates + 2;              % increase number of states by 2 (1 for xi, 1 for yi)
                [xest,Pest] = updateNew(zlaser(:,i));   % update state and covariance
            else                                        % not 1st landmark -> check to make sure no LMs in range
                [dist,index]=nearby_LMs(zlaser(:,i));   % dist from observation to already incorporated LMs (< 3m)
                min_dist = 2;                           % minimum distance between landmark and observation
                if dist>min_dist                        % if dist to nearby_LM is > "min_dist", add observation as new LM
                    new_state(zlaser(:,i));             % add new state (i.e. increase all matrix sizes)
                    numStates = numStates + 2;          % increase number of states by 2 (1 for xi, 1 for yi)
                    [xest,Pest] = updateNew(zlaser(:,i));               % update state and covariance
                else                                                    % else find closest LM and associate
                    closest = 4 + 2*(index-1);                          % find LM which is closest to observation
                    [xest,Pest] = updateExisting(zlaser(:,i),closest);  % update state and covariance
                end
            end
        end
    end
    % LB_add: display information
    if( floor(t/tfinal*100) > n_dsp)
        clc;
        n_dsp = floor(t/tfinal*100);
        str = sprintf('运行卡尔曼滤波... [%d%c]', n_dsp,'%');
        disp(str);
    end
end  % end while
%-------------- end of  EXTENDED KALMAN FILTER PROCEDURE -----------------
 
%--------------------------- PLOT RESULTS --------------------------------
figure(1); hold on;
plot(GPSLon(1:560),GPSLat(1:560));          % robot position (true)
plot(xest(1,1:k),xest(2,1:k),'g');          % robot position (model)
plot(landmarks(:,1),landmarks(:,2),'*');    % landmark position (true)
for i=4:2:numStates
    plot(xest(i,k),xest(i+1,k),'g*');     % landmark position (model)
end
legend('机器人实际位置','机器人slam位置','路标实际位置','路标slam位置');
xlabel('x [meters]'); ylabel('y [meters]');
axis([-10 20 -25 20]);
 
%% Estimate error(估计误差)  
x_error1 = GPSLon(1:560)-xest(1,1:560);
x_error2 = GPSLat(1:560)-xest(2,1:560);
for i=4:2:numStates 
    landmarks_error1=landmarks(:,2)-xest(i+1,k);
end
for i=4:2:numStates 
    landmarks_error2=landmarks(:,1)-xest(i,k);
end
%% Graph 2  
figure(2)   
plot(x_error1),grid on;  
title('机器人位置误差 on X axis')  
xlabel('时间 [sec]')  
ylabel('位置均方根误差 [m]') 
axis([0 112 -25 25]);
hold off; 
%% Graph 3  
figure(3)   
plot(x_error2),grid on;  
title('机器人位置误差 on Y axis')  
xlabel('时间 [sec]')  
ylabel('位置均方根误差 [m]') 
axis([0 112 -25 25]);
hold off;
%% Graph 4  
figure(4)   
plot(landmarks_error1),grid on;  
title('路标位置误差 on x axis')  
xlabel('时间 [sec]')  
ylabel('位置均方根误差 [m]') 
axis([0 18 -25 25]);
hold off;
%% Graph 5 
figure(5)   
plot(landmarks_error2),grid on;  
title('路标位置误差 on Y axis')  
xlabel('时间 [sec]')  
ylabel('位置均方根误差 [m]') 
axis([0 18 -25 25]);
hold off;
%------------------------ end of PLOT RESULTS -----------------------------
