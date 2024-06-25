close all; clear; clc;
addpath ./tasks;

% Environment Configurations
tx_node_number = 1;          % Number of Tx users
rx_node_number = 2;          % Number of Rx users
analog_antenna_number = 16;  % Number of Tx antennas of analog
digital_antenna_number = 2;  % Number of Tx antennas of digital
rx_antenna_number = 1;       % Number of Rx antennas
codebook_size = 10;
% Generate receivers with beam direction and distances
% [x y] = rand_rx_location_list(tx_beam_direction_idx)
rand_rx_location_list = [];
for tx_beam = 0:10:180
    tmp = [];
    for d = 50:50:500
        offset = -5 + 10 * rand();       % -5~5 degrees
        x = d * cosd(tx_beam + offset);  % Add a small random offset
        y = d * sind(tx_beam + offset);  % Add a small random offset
        tmp = [tmp x y];
    end
    rand_rx_location_list = [rand_rx_location_list; tmp];
end

% Define coordination and power for transmitter
origin = [0, 0];
tx_location = origin;
P_tx_dBm = 10;          % Transmission power of Tx (dBm)
N0_dBm = -95;           % Assume noise power is -90 dBm

% Randomly select receiver's coordination at d=200
rand_idx = randperm(numel(0:10:180), 2);
rx1_location = rand_rx_location_list(rand_idx(1), (7:8));
rx2_location = rand_rx_location_list(rand_idx(2), (7:8));
% disp(rand_idx)
%disp(rand_rx_location_list)
%disp(rx1_location)
%disp(rx2_location)
%disp(rand_idx)

% TODO: Implement Analog Beamforming functions in /tasks
% Hint: you can adjust input/output for reports
[rx1_SNR_dbm, rx2_INR_dbm] = analog_beamforming(codebook_size, P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, analog_antenna_number);

fprintf('Task1: SNR Calculation\n');
fprintf('Receiver1\tSNR: %f dBm\tReceiver2 INR: %f dBm\n', rx1_SNR_dbm, rx2_INR_dbm);

% calculate average SNR and INR through distance [50:50:500]

% the array to store the average value of SNR and INR
avg_snr_8 = zeros(10);
avg_inr_8 = zeros(10);
avg_snr_16 = zeros(10);
avg_inr_16 = zeros(10);

for antenna_amount = 8:8:16
for dist = 50:50:500
    analog_antenna_number = antenna_amount;
    snr_sum = 0;
    inr_sum = 0;
    times = 0;

    for repeat = 1:1:10
    rand_rx_location_list1 = [];
    for tx_beam = 0:10:180
        tmp = [];
        for d = 50:50:500
            offset = -5 + 10 * rand();       % -5~5 degrees
            x = d * cosd(tx_beam + offset);  % Add a small random offset
            y = d * sind(tx_beam + offset);  % Add a small random offset
            tmp = [tmp x y];
        end
        rand_rx_location_list1 = [rand_rx_location_list1; tmp];
    end

    % generate the random location when distance is dist
    rand_idx = randperm(numel(0:10:180), 2);
    rx1_location = rand_rx_location_list(rand_idx(1), ((dist*2)/50-1:(dist*2)/50));
    rx2_location = rand_rx_location_list(rand_idx(2), ((dist*2)/50-1:(dist*2)/50));

    % use analog beamforming func. to get SNR of user 1 and INR of user 2
    [pathloss, rx1_SNR_dbm, rx2_INR_dbm] = analog_beamforming(codebook_size, P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, analog_antenna_number);
    
    % sum up the snr and inr
    snr_sum = snr_sum + rx1_SNR_dbm;
    inr_sum = inr_sum + rx2_INR_dbm;
    %timmes = times + 1;
    end

    avg_snr = snr_sum/10;
    avg_inr = inr_sum/10;
    if antenna_amount == 8
        avg_snr_8((dist/50)) = avg_snr;
        avg_inr_8((dist/50)) = avg_inr;
    else
        avg_snr_16((dist/50)) = avg_snr;
        avg_inr_16((dist/50)) = avg_inr;
    end

    fprintf('Antenna Numbers: %f,Distance: %f m, average SNR: %f dBm, average INR: %f dBm\n', antenna_amount, dist, avg_snr, avg_inr);
end
end

% plot the picture (# of antenna: 8)
figure
x = linspace(50, 500, 10);
y = avg_snr_8((x/50));
plot(x, y)

xlabel('distance(m)')
ylabel('dBm')

hold on
y2 = avg_inr_8((x/50));
plot(x, y2)
legend('SNR', 'INR')
hold off
% antenna number = 16
figure
x = linspace(50, 500, 10);
y = avg_snr_16((x/50));
plot(x, y)

xlabel('distance(m)')
ylabel('dBm')

hold on
y2 = avg_inr_16((x/50));
plot(x, y2)
legend('SNR', 'INR')
hold off


% output the result of SNR and INR of 10 topologies when d=200m, # of
% antenna = 16
snr_task1_2 = zeros(10);
inr_task1_2 = zeros(10);
for repeat = 1:1:10
    rand_rx_location_list = [];
    for tx_beam = 0:10:180
        tmp = [];
        for d = 50:50:500
            offset = -5 + 10 * rand();       % -5~5 degrees
            x = d * cosd(tx_beam + offset);  % Add a small random offset
            y = d * sind(tx_beam + offset);  % Add a small random offset
            tmp = [tmp x y];
        end
        rand_rx_location_list = [rand_rx_location_list; tmp];
    end

    % Randomly select receiver's coordination at d=200
    rand_idx = randperm(numel(0:10:180), 2);
    rx1_location = rand_rx_location_list(rand_idx(1), (7:8));
    rx2_location = rand_rx_location_list(rand_idx(2), (7:8));
    

    analog_antenna_number = 16;
    [pathloss, rx1_SNR_dbm, rx2_INR_dbm] = analog_beamforming(codebook_size, P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, analog_antenna_number);
    snr_task1_2(repeat) = rx1_SNR_dbm;
    inr_task1_2(repeat) = rx2_INR_dbm;

    fprintf('Topology: %d, SNR: %f dBm, INR: %f dBm\n', repeat, snr_task1_2(repeat), inr_task1_2(repeat))
end

figure
x = linspace(1, 10, 10);
y = snr_task1_2(x);
plot(x, y)

xlabel('times')
ylabel('dBm')
title('Task 1-2')

hold on
y2 = inr_task1_2(x);
plot(x, y2)
legend('SNR', 'INR')
hold off


% plot the Prx,1 (in dBm) of 10 topologies for various codebook sizes (19,
% 37, 73, i.e. [0:10:180], [0:5:180], [0:2.5:180]) when d=200m, antenna
% number = 16

% [0:?:180]
% you can adjust ? by changing the tx_beam_direction in
% analog_beamforming.m
ptx_dbm10 = zeros(10);
ptx_dbm5 = zeros(10);
ptx_dbm25 = zeros(10);
for repeat = 1:1:10
    rand_rx_location_list = [];
    for tx_beam = 0:10:180
        tmp = [];
        for d = 50:50:500
            offset = -5 + 10 * rand();       % -5~5 degrees
            x = d * cosd(tx_beam + offset);  % Add a small random offset
            y = d * sind(tx_beam + offset);  % Add a small random offset
            tmp = [tmp x y];
        end
        rand_rx_location_list = [rand_rx_location_list; tmp];
    end

    % Randomly select receiver's coordination at d=200
    rand_idx = randperm(numel(0:10:180), 2);
    rx1_location = rand_rx_location_list(rand_idx(1), (7:8));
    rx2_location = rand_rx_location_list(rand_idx(2), (7:8));
    

    analog_antenna_number = 16;
    codebook_size = 10;
    [pathloss_rx_dbm, rx1_SNR_dbm, rx2_INR_dbm] = analog_beamforming(codebook_size, P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, analog_antenna_number);
    ptx_dbm10(repeat) = pathloss_rx_dbm;

    %fprintf('Codebook size: [0:10:180], Prx,1: %f dBm\n', ptx_dbm10(repeat))

    codebook_size = 5;
    [pathloss_rx_dbm, rx1_SNR_dbm, rx2_INR_dbm] = analog_beamforming(codebook_size, P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, analog_antenna_number);
    ptx_dbm5(repeat) = pathloss_rx_dbm;

    %fprintf('Codebook size: [0:5:180], Prx,1: %f dBm\n', ptx_dbm5(repeat))

    codebook_size = 2.5;
    [pathloss_rx_dbm, rx1_SNR_dbm, rx2_INR_dbm] = analog_beamforming(codebook_size, P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, analog_antenna_number);
    ptx_dbm25(repeat) = pathloss_rx_dbm;

    %fprintf('Codebook size: [0:2.5:180], Prx,1: %f dBm\n', ptx_dbm25(repeat))
end

for repeat = 1:1:10
    fprintf('Codebook size: [0:10:180], Prx,1: %f dBm\n', ptx_dbm10(repeat))
end

for repeat = 1:1:10
    fprintf('Codebook size: [0:5:180], Prx,1: %f dBm\n', ptx_dbm5(repeat))
end

for repeat = 1:1:10
    fprintf('Codebook size: [0:2.5:180], Prx,1: %f dBm\n', ptx_dbm25(repeat))
end

figure
x = linspace(1, 10, 10);
y = ptx_dbm10(x);
plot(x, y)

xlabel('times')
ylabel('dBm')
title('Codebook')
legend('Prx,1')

hold on
y2 = ptx_dbm5(x);
plot(x, y2)
y3 = ptx_dbm25(x);
plot(x, y3)
legend('[0:10:180]', '[0:5:180]', '[0:2.5:180]')
hold off


% TODO: Implement Digital Beamforming functions in /tasks
% Hint: you can adjust input/output for reports
[h_eq, error_rx1_w, rx1_SNR_dbm, rx2_SNR_dbm, ori_rx1_SNR_dbm, ori_rx2_SNR_dbm] = digital_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, tx_node_number*digital_antenna_number, rx_node_number*rx_antenna_number);

fprintf('Task2: SNR Calculation\n');
fprintf('Receiver1\tSNR: %f dBm\t without precoding SNR: %f dBm\n', rx1_SNR_dbm, ori_rx1_SNR_dbm);
fprintf('Receiver2\tSNR: %f dBm\t without precoding SNR: %f dBm\n', rx2_SNR_dbm, ori_rx2_SNR_dbm);

% plot the average SNR of two users w/ and w/o ZFBF (x-axis: distances)
% d = [50:50:500]
avg_rx1_snr_w = zeros(10);
avg_rx2_snr_w = zeros(10);
avg_rx1_snr_wo = zeros(10);
avg_rx2_snr_wo = zeros(10);

for dist = 50:50:500
    sum_rx1_snr_w = 0;
    sum_rx2_snr_w = 0;
    sum_rx1_snr_wo = 0;
    sum_rx2_snr_wo = 0;
    for i = 1:1:10
        rand_rx_location_list1 = [];
        for tx_beam = 0:10:180
            tmp = [];
            for d = 50:50:500
                offset = -5 + 10 * rand();       % -5~5 degrees
                x = d * cosd(tx_beam + offset);  % Add a small random offset
                y = d * sind(tx_beam + offset);  % Add a small random offset
                tmp = [tmp x y];
            end
            rand_rx_location_list1 = [rand_rx_location_list1; tmp];
        end

        % generate the random location when distance is dist
        rand_idx = randperm(numel(0:10:180), 2);
        rx1_location = rand_rx_location_list(rand_idx(1), ((dist*2)/50-1:(dist*2)/50));
        rx2_location = rand_rx_location_list(rand_idx(2), ((dist*2)/50-1:(dist*2)/50));
        [h_eq, error_rx1_w, rx1_SNR_dbm, rx2_SNR_dbm, ori_rx1_SNR_dbm, ori_rx2_SNR_dbm] = digital_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, tx_node_number*digital_antenna_number, rx_node_number*rx_antenna_number);

        sum_rx1_snr_w = sum_rx1_snr_w + rx1_SNR_dbm;
        sum_rx2_snr_w = sum_rx2_snr_w + rx2_SNR_dbm;
        sum_rx1_snr_wo = sum_rx1_snr_wo + ori_rx1_SNR_dbm;
        sum_rx2_snr_wo = sum_rx2_snr_wo + ori_rx2_SNR_dbm;
    end

    avg_rx1_snr_w((dist/50)) = sum_rx1_snr_w/10;
    avg_rx2_snr_w((dist/50)) = sum_rx2_snr_w/10;
    avg_rx1_snr_wo((dist/50)) = sum_rx1_snr_wo/10;
    avg_rx2_snr_wo((dist/50)) = sum_rx2_snr_wo/10;
    fprintf('Distance: %f m, SNR1 w/: %f dBm, SNR2 w/: %f dBm, SNR1 w/o: %f dBm, SNR2 w/o: %f dBm\n', dist, avg_rx1_snr_w((dist/50)), avg_rx2_snr_w((dist/50)), avg_rx1_snr_wo((dist/50)), avg_rx2_snr_wo((dist/50)));
end

figure
x = linspace(50, 500, 10);
y = avg_rx1_snr_w((x/50));
plot(x, y);

xlabel('distance (m)')
ylabel('dBm')
title('Task 2-1 (user 1)')

hold on
y2 = avg_rx1_snr_wo((x/50));
plot(x, y2)
legend('w/ ZFBF', 'w/o ZFBF')
hold off

figure
x = linspace(50, 500, 10);
y = avg_rx2_snr_w((x/50));
plot(x, y);

xlabel('distance (m)')
ylabel('dBm')
title('Task 2-1 (user 2)')

hold on
y2 = avg_rx2_snr_wo((x/50));
plot(x, y2)
legend('w/ ZFBF', 'w/o ZFBF')
hold off

% plot h_eq and error (in dBm) of R1 with ZFBF when d=200m
rand_rx_location_list2 = [];
for tx_beam = 0:10:180
    tmp = [];
    for d = 50:50:500
        offset = -5 + 10 * rand();       % -5~5 degrees
        x = d * cosd(tx_beam + offset);  % Add a small random offset
        y = d * sind(tx_beam + offset);  % Add a small random offset
        tmp = [tmp x y];
    end
    rand_rx_location_list2 = [rand_rx_location_list2; tmp];
end
% Randomly select receiver's coordination at d=200
rand_idx = randperm(numel(0:10:180), 2);
rx1_location = rand_rx_location_list2(rand_idx(1), (7:8));
rx2_location = rand_rx_location_list2(rand_idx(2), (7:8));

h_eq_rx1 = zeros([1, 10]);
error_rx1_with = zeros([1, 10]);

for i = 1:1:10
    [h_eq, error_rx1_w, rx1_SNR_dbm, rx2_SNR_dbm, ori_rx1_SNR_dbm, ori_rx2_SNR_dbm] = digital_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, tx_node_number*digital_antenna_number, rx_node_number*rx_antenna_number);
    h_eq_rx1(1, i) = h_eq;
    error_rx1_with(1, i) = error_rx1_w;
end
%disp(h_eq_rx1)

%disp(error_rx1_with)
figure
x = linspace(1, 10, 10);
y = h_eq_rx1(x);
plot(x, y)

figure
x = linspace(1, 10, 10);
y = error_rx1_with(x);
plot(x, y)




