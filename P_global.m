function [x_k_average_P, y_k_average_P, average_A_P] = P_global(n , p , t_end , x_0 , n_selfish , ref , Kp)

%Function giving the following output of a P-controlled system with global
%visibility:
% - plot of opinion dynamics in the P-controlled system
% - opinion values evolution
% - network average evolution

% Reference sequence (for plots)
ref_seq = mean(x_0(n_selfish+1:end)) * ones(t_end+1 , 1);

% Computing the average Mask
average_mask = p/(n-1)*(ones(n)-eye(n))+eye(n);

% State-Space Representation 
average_A = (1/(p+1))*average_mask;
B = [eye(n_selfish) ; zeros(n - n_selfish , n_selfish)];
C = (1/n)*[ones(n_selfish , n)];

% Open-Loop System
average_A_P = average_A - B*Kp*C;

% Definition of variables for storage of state evolution and Measurements
% evolution
x_k_average_P = zeros(n , t_end+1); 
y_k_average_P = zeros(n_selfish , t_end+1); 

x_k_average_P(: , 1) = x_0;
y_k_average_P(: , 1) = C * x_0;

% Building the sequence of state matrices for P-Controlled Closed Loop 
for k  = 1:t_end
    x_k_average_P(: , k+1) = average_A_P * x_k_average_P(: , k) + B*Kp*ref;
    y_k_average_P(: , k+1) = C * x_k_average_P(: , k+1);
end

%Plotting Opinion Dynamics
figure(100) ;  hold on;
plot(0:1:t_end , x_k_average_P(1:3 ,:) ,  'LineWidth' , 1.5); hold on;
plot(0:1:t_end, x_k_average_P(n_selfish+1 ,:),  'LineWidth' , 1.5); 
plot(0:1:t_end , y_k_average_P(1 ,:) , 'LineWidth' , 1.5);
plot(0:1:t_end, ref_seq, 'k .' , 'MarkerSize' , 1.2);
legend( 'Coordinator 1' ,'Coordinator 2' ,'Coordinator 3' , 'Standard Agent 1' , 'Network average' , 'Reference'  ,'Location' ,'SouthEast');
title('Global, NO saturation, Mean reference, P');
hold off;

end

