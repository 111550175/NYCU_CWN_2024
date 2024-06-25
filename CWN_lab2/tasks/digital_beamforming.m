function [h_eq, return_avg_error_rx1_w_dbm, rx1_SNR_dbm, rx2_SNR_dbm, rx1_noW_SNR_dbm, rx2_noW_SNR_dbm] = digital_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, tx_nums, rx_nums)
addpath ./ewa_function;

% Calculate P_rx
% Hint: Use friis equation
freq = 2.4e9;
PL1_dBm = friis_equation(freq, 1, 1, norm(rx1_location - tx_location));
PL2_dBm = friis_equation(freq, 1, 1, norm(rx2_location - tx_location));
P_rx1_dBm = P_tx_dBm + PL1_dBm;
P_rx2_dBm = P_tx_dBm + PL2_dBm;
%disp(P_rx1_dBm)
%disp(P_rx2_dBm)

P_rx1_mWatt = 10^((P_rx1_dBm)/10);
P_rx2_mWatt = 10^((P_rx2_dBm)/10);

% Generate random channel h~N(0,1)
H = (randn(rx_nums, tx_nums) + randn(rx_nums, tx_nums) * 1i) ./ sqrt(2) .* sqrt([P_rx1_mWatt; P_rx2_mWatt]);

% TODO1: Generate precoded weight W (h')
% Hint1: You can reference the equation for calculating W on page 35 of PPT L9
% Hint2: Or you can just use the MATLAB api to get inverse matrix of H
W = ones(size(H));
%disp(size(H))
%disp(H)
% use the api of matlab to calculate the inverse of H
W = inv(H);

% TODO2: Scaled W into unit power (power = 1)
% sum up all the element in W, then take a square root to scaled it
% correctly
scaled_w = sqrt(sum(abs(W).^2, "all"));
%disp(scaled_w)
W = (1/(scaled_w)) * W;
%disp(W)

% Hint: Uncomment the equation below to verify if W is converted to unit power (W_power = 1)
W_power = sum(abs(W).^2, 'all');
%fprintf('W_power: %f\n', W_power);

% Generate random transmitted signals
num_data = 1000;
x = randn(rx_nums, num_data);
%disp(x)

% Generate random noises
N0_mWatt = 10^((N0_dBm)/10); % Convert noise power to mW
n = (randn(rx_nums, num_data) + randn(rx_nums, num_data) * 1i) ./ sqrt(2) * sqrt(N0_mWatt);

% Generate receive signals 
y_noW = H * x + n; % without precoding
y = H * W * x + n; % with precoding

% TODO3: Calculate H_eq with & without ZFBF
% Hint: without ZFBF: summmation one rx channel by all tx channel. H_eq = sum H_rx,i i=tx
% Hint: with ZFBF: H * W = I * constant. Please calculate the constant (H_eq)
% H_eq without ZFBF:
heq_wo_1 = H(1, 1) + H(1, 2);
heq_wo_2 = H(2, 1) + H(2, 2);

% H_eq with ZFBF: (you can get the constant by the index(1, 1) of the
% multiplication of H and W)
I = H * W;
%disp(I)
h_eq = I(1, 1);
%disp(h_eq)

% TODO4: Decode signal with both with & without ZFBF H_eq
% Hint: Estimate x_hat based on the received signal y
%disp(y_noW)
% calculate x-hat by y/h_eq (w/o ZFBF)
x_hat_rx1_wo = zeros([1, 1000]);
x_hat_rx2_wo = zeros([1, 1000]);
for i = 1:1:1000
    x_hat_rx1_wo(1, i) = y_noW(1, i)/heq_wo_1;
    x_hat_rx2_wo(1, i) = y_noW(2, i)/heq_wo_2;
end

% calculate x-hat by y/h_eq (w/ ZFBF)
x_hat_rx1_w = zeros([1, 1000]);
x_hat_rx2_w = zeros([1, 1000]);
for i = 1:1:1000
    x_hat_rx1_w(1, i) = y(1, i)/h_eq;
    x_hat_rx2_w(1, i) = y(2, i)/h_eq;
end

%disp(x_hat_rx2_w)

% TODO5: Calculate decoding errors (|x-x'|^2) and SNR (|x|^2/|x-x'|^2)
% Hint: The difference between the transmitted and estimated signals
rx1_SNR_dbm = 0;
rx2_SNR_dbm = 0;
rx1_noW_SNR_dbm = 0;
rx2_noW_SNR_dbm = 0;

% calculate decoding error (w/o):
% rx1
error_rx1 = 0;
for i = 1:1:1000
    %error_rx1 = error_rx1 + (abs(x(1, i) - x_hat_rx1_wo(1, i)))^2;
    error_rx1 = error_rx1 + ( real( x(1, i) - x_hat_rx1_wo(1, i) ) )^2 + ( imag(x(1, i) - x_hat_rx1_wo(1, i)) )^2;
end

avg_error_rx1 = error_rx1/1000;
%rx2
error_rx2 = 0;
for i = 1:1:1000
    %error_rx2 = error_rx2 + (abs(x(2, i) - x_hat_rx2_wo(1, i)))^2;
    error_rx2 = error_rx2 + ( real( x(2, i) - x_hat_rx2_wo(1, i) ) )^2 + ( imag(x(2, i) - x_hat_rx2_wo(1, i)) )^2;
end

avg_error_rx2 = error_rx2/1000;

% calculate SNR (w/o):
% rx1
sound_rx1 = 0;
for i = 1:1:1000
    %sound_rx1 = sound_rx1 + (abs(x(1, i)))^2;
    sound_rx1 = sound_rx1 + ( real(x(1, i)) )^2 + ( imag(x(1, i)) )^2;
end
avg_sound_rx1 = sound_rx1/1000;
%rx2
sound_rx2 = 0;
for i = 1:1:1000
    %sound_rx2 = sound_rx2 + (abs(x(2, i)))^2;
    sound_rx2 = sound_rx2 + ( real(x(2, i)) )^2 + ( imag(x(2, i)) )^2;
end
avg_sound_rx2 = sound_rx2/1000;
% first take the mean of sound and noise, then divide them to get avg SNR
rx1_noW_SNR_dbm = 10*log10(avg_sound_rx1/avg_error_rx1);
rx2_noW_SNR_dbm = 10*log10(avg_sound_rx2/avg_error_rx2);

% w/ ZFBF:
% calculate decoding error:
% rx1
error_rx1_w = 0;
for i = 1:1:1000
    %error_rx1_w = error_rx1_w + (abs(x(1, i) - x_hat_rx1_w(1, i)))^2;
    error_rx1_w = error_rx1_w + ( real( x(1, i) - x_hat_rx1_w(1, i) ) )^2 + ( imag(x(1, i) - x_hat_rx1_w(1, i)) )^2;
end
avg_error_rx1_w = error_rx1_w/1000;

%rx2
error_rx2_w = 0;
for i = 1:1:1000
    %error_rx2_w = error_rx2_w + (abs(x(2, i) - x_hat_rx2_w(1, i)))^2;
    error_rx2_w = error_rx2_w + ( real( x(2, i) - x_hat_rx2_w(1, i) ) )^2 + ( imag(x(2, i) - x_hat_rx2_w(1, i)) )^2;
end
avg_error_rx2_w = error_rx2_w/1000;

return_avg_error_rx1_w_dbm = 10*log10(avg_error_rx1_w/0.001);

%calculate SNR:
rx1_SNR_dbm = 10*log10(avg_sound_rx1/avg_error_rx1_w);
rx2_SNR_dbm = 10*log10(avg_sound_rx2/avg_error_rx2_w);

%disp(avg_sound_rx1)
%disp(avg_sound_rx2)
%disp(avg_error_rx1_w)
%disp(avg_error_rx2_w)
%disp(avg_error_rx1)
%disp(avg_error_rx2)

% handle the warning
h_eq = real(h_eq);
