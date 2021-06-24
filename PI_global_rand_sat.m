function [x_k_PI , y_k_PI] = PI_global_rand_sat(n , p , t_end , x_0 , n_selfish , ref , A_sequence , Kp , Ki);


% We are checking if the controller works when using the random sequence of
% i.i.d. matrices intead of their expected value, still with global
% visibility

%Function giving the following output of a PI-controlled system with global
%visibility:
% - plot of opinion dynamics in the PI-controlled system
% - opinion values evolution
% - network average evolution

% Reference sequence (for plots)
ref_seq = mean(x_0(n_selfish+1:end)) * ones(t_end+1 , 1);

% State-Space Representation 
A_PI_sequence = zeros(n+n_selfish , n+n_selfish , t_end);
B = [eye(n_selfish) ; zeros(n - n_selfish , n_selfish)];
C = (1/n)*[ones(n_selfish , n)];

% Definition of variables for storage of the augmented state evolution and the Measurements
% evolution
x_0_aug=[x_0; zeros(n_selfish,1)];
x_k_PI = zeros(n + n_selfish , t_end+1);
y_k_PI = zeros(n_selfish , t_end+1);
x_k_PI(1:n , 1) = sat_function(x_0_aug(1:n , 1));
y_k_PI(: , 1) = C*x_0;

% Building the sequence of state matrices for PI-Controlled Closed Loop 
for k  = 1:t_end
    A_PI_sequence(: , : , k) = [A_sequence(: , : , k)-B*Kp*C , B*Ki ; -C , eye(n_selfish)];
    x_k_PI(1:n , k+1) = sat_function( A_PI_sequence(1:n , : , k) * x_k_PI(: , k) + B*Kp*ref );
    x_k_PI(n+1:end , k+1) = -C * x_k_PI(1:n , k) + eye(n_selfish)*x_k_PI(n+1:end , k);
    y_k_PI(: , k+1) = C * x_k_PI(1:n , k+1);
end

%Plotting opinion Dynamics
figure(2000) ;  hold on;
plot(0:1:t_end , x_k_PI(1:3 ,:) ,  'LineWidth' , 1.5); hold on;
plot(0:1:t_end, x_k_PI(n_selfish+1 ,:),  'LineWidth' , 1.5);
% plot(0:1:t_end, x_k_PI(n+1 ,:),  'LineWidth' , 1.5);
plot(0:1:t_end , y_k_PI(1 ,:) , 'LineWidth' , 1.5);
plot(0:1:t_end, ref_seq, 'k .' , 'MarkerSize' , 1.1);
legend( 'Coordinator 1' ,'Coordinator 2' ,'Coordinator 3' , 'Standard Agent 1' , 'Integrated Error Coordinator 1' , 'Network average' , 'Reference');
title('Global, saturation, Mean reference, PI , random sequence');
hold off;


end

