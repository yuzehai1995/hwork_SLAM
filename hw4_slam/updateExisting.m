% This is the measurement update step of the Extended Kalman Filter.  It is used
% to update an *existing* landmark and the state and covariance matrices.
function [xest, Pest] = updateExisting(zlaser,landmark)

    global xest;        % state matrix
    global Pest;        % covariance matrix
    global numStates;   % number of states before entering this loop
    global k;    
    
    dx = xest(landmark,k) - xest(1,k);      % distance between existing landmark and x-pos of vehicle
    dy = xest(landmark+1,k) - xest(2,k);    % distance between existing landmark and y-pos of vehicle
    r = sqrt(dx^2 + dy^2);
    
    % Sensor Model - range and bearing
    H = [r; atan2(dy,dx) - xest(3,k) + pi/2];
    H(2) = normalizeAngle(H(2));
  
    % Calculate Jacobian Jh (dh/dx)
    Jh(1,1) = -dx/r;   Jh(1,2) = -dy/r;   Jh(1,3) = 0;
    Jh(2,1) = dy/(r*r);   Jh(2,2) = -dx/(r*r);   Jh(2,3) = -1;
    if(numStates>3)
        Jh(1:2,4:numStates) = zeros(2,numStates-3);
    end
    Jh(1,landmark) = dx/r;   Jh(1,landmark+1) = dy/r;
    Jh(2,landmark) = -dy/(r*r);     Jh(2,landmark+1) = dx/(r*r);
    % --LB_debug: changed from Jh(1,..)=-dx/r, four operator are all reversed
    % --LB_debug: this is the biggest error in the debug procedure
  
    V = eye(2);
    R = [0.01, 0; 0, (pi/180)^2];
    
    innov = zlaser - H;
    innov(2) = normalizeAngle(innov(2));
               
    K = Pest * Jh' / (Jh*Pest*Jh'+V*R*V');             % Kalman gain
    xest(:,k) = xest(:,k) + K*innov;     % update state matrix
    xest(3,k) = normalizeAngle(xest(3,k));      % keep phi between -180 and 180
    I = eye(size(K*Jh));
    Pest = (I-K*Jh)*Pest;          % udpate covariance matrix
   
return;