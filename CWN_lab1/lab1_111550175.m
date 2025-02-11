%task1 - Prx
%testing
count =1;
y1 = zeros(12, 1);
y2 = zeros(12, 1);
y3 = zeros(12, 1);
y4 = zeros(12, 1);

rng(0)
a = randn;
b = randn;
ranBN = randi([0, 1], 1, 300000);
ori_awgn = zeros(300000, 1);

bpsk_xi_prime = zeros(300000, 1);
qpsk_xi_prime = zeros(150000, 1);
qam_16_xi_prime = zeros(75000, 1);
qam_64_xi_prime = zeros(50000, 1);

for i = 1:300000
    ori_noise_x = randn();
    ori_noise_y = randn();
    ori_awgn(i) = complex(ori_noise_x, ori_noise_y);
    %scaled_noise = (n0) * complex(ori_noise_x, ori_noise_y);
    %awgn(i) = scaled_noise;
end

for distance = 50:50:600
t1 = ['Distance: ', num2str(distance)];
%y1 = zeros(12, 1);
%y2 = zeros(12, 1);
%y3 = zeros(12, 1);
%y4 = zeros(12, 1);

disp(t1)
tx_gain = 1;
rx_gain = 1;
tx_power_dbm = 10;
wavelength = 3e+8/2.4e+9;

rx_power_dbm = tx_power_dbm + 10*log10(tx_gain) + 10*log10(rx_gain) + 20*log10(wavelength/(4*pi*distance));
x = ['Prx in dBm: ', num2str(rx_power_dbm)];
disp(x)

rx_power_watt = (10^(rx_power_dbm/10))*(1e-3);
y = ['Prx in Watt: ', num2str(rx_power_watt)];
disp(y)

% task2 - channel generation
%a = randn;
%b = randn;

h = complex(a, b);
scaling_h = h * ((rx_power_watt)^(1/2));
% disp(num2str(h))

% task3 - modulation

% generate 300,000 random bits
%ranBN = randi([0, 1], 1, 300000);

%BPSK
bpsk = zeros(300000, 1);

const_point_bpsk = [complex(1, 0), complex(-1, 0)];
bpsk_comp = zeros(300000, 1);

for i = 1:300000
    bpsk(i) = ranBN(i);
    bpsk_comp(i) = const_point_bpsk(bpsk(i)+1);
end

%for i = 1:10
    %disp((bpsk_comp(i)))
%end

%QPSK
qpsk = zeros(150000, 1);
const_point_qpsk = [complex(1/sqrt(2), 1/sqrt(2)), complex(1/sqrt(2), -1/sqrt(2)),...
    complex(-1/sqrt(2), 1/sqrt(2)), complex(-1/sqrt(2), -1/sqrt(2))];
qpsk_comp = zeros(150000, 1);
for i = 1:150000
    sum = 0;
    for j = 1:2
        sum = sum + ranBN((i-1)*2+(j-1)+1);
        sum = sum * 2;
    end
    sum = sum/2;
    qpsk(i) = sum;
    qpsk_comp(i) = const_point_qpsk(qpsk(i)+1);
    % disp(num2str(ranBN(i)))
end
%for i = 1:10
    %disp(qpsk_comp(i))
%end

%16-QAM
qam_16 = zeros(75000, 1);
const_point_qam_16 = [complex(3/sqrt(10),3/sqrt(10)), complex(3/sqrt(10),1/sqrt(10)), complex(3/sqrt(10), -1/sqrt(10)),... 
    complex(3/sqrt(10),-3/sqrt(10)), complex(1/sqrt(10),3/sqrt(10)), complex(1/sqrt(10),1/sqrt(10)),...
    complex(1/sqrt(10),-1/sqrt(10)), complex(1/sqrt(10),-3/sqrt(10)), complex(-1/sqrt(10),3/sqrt(10)),...
    complex(-1/sqrt(10),1/sqrt(10)), complex(-1/sqrt(10),-1/sqrt(10)), complex(-1/sqrt(10),-3/sqrt(10)),...
    complex(-3/sqrt(10),3/sqrt(10)), complex(-3/sqrt(10),1/sqrt(10)), complex(-3/sqrt(10),-1/sqrt(10)),...
    complex(-3/sqrt(10),-3/sqrt(10))];

qam_16_comp = zeros(75000, 1);
for i = 1:75000
    sum = 0;
    for j = 1:4
        sum = sum + ranBN((i-1)*4+(j-1)+1);
        sum = sum*2;
    end
    sum = sum/2;
    qam_16(i) = sum;
    qam_16_comp(i) = const_point_qam_16(qam_16(i)+1);
end
%for i = 1:10
    %disp(num2str(qam_16_comp(i)))
%end

%64-QAM
qam_64 = zeros(50000, 1);
const_point_qam_64 = [complex(7/sqrt(42),7/sqrt(42)), complex(7/sqrt(42),5/sqrt(42)), complex(7/sqrt(42),3/sqrt(42)),...
    complex(7/sqrt(42),1/sqrt(42)), complex(7/sqrt(42),-1/sqrt(42)), complex(7/sqrt(42),-3/sqrt(42)),...
    complex(7/sqrt(42),-5/sqrt(42)), complex(7/sqrt(42),-7/sqrt(42)), complex(5/sqrt(42),7/sqrt(42)),...
    complex(5/sqrt(42),5/sqrt(42)), complex(5/sqrt(42),3/sqrt(42)), complex(5/sqrt(42),1/sqrt(42)),...
    complex(5/sqrt(42),-1/sqrt(42)), complex(5/sqrt(42),-3/sqrt(42)), complex(5/sqrt(42),-5/sqrt(42)),...
    complex(5/sqrt(42),-7/sqrt(42)), complex(3/sqrt(42),7/sqrt(42)), complex(3/sqrt(42),5/sqrt(42)),...
    complex(3/sqrt(42),3/sqrt(42)), complex(3/sqrt(42),1/sqrt(42)), complex(3/sqrt(42),-1/sqrt(42)),...
    complex(3/sqrt(42),-3/sqrt(42)), complex(3/sqrt(42),-5/sqrt(42)), complex(3/sqrt(42),-7/sqrt(42)),...
    complex(1/sqrt(42),7/sqrt(42)), complex(1/sqrt(42),5/sqrt(42)), complex(1/sqrt(42),3/sqrt(42)),...
    complex(1/sqrt(42),1/sqrt(42)), complex(1/sqrt(42),-1/sqrt(42)), complex(1/sqrt(42),-3/sqrt(42)),...
    complex(1/sqrt(42),-5/sqrt(42)), complex(1/sqrt(42),-7/sqrt(42)), complex(-1/sqrt(42),7/sqrt(42)),...
    complex(-1/sqrt(42),5/sqrt(42)), complex(-1/sqrt(42),3/sqrt(42)), complex(-1/sqrt(42),1/sqrt(42)),...
    complex(-1/sqrt(42),-1/sqrt(42)), complex(-1/sqrt(42),-3/sqrt(42)), complex(-1/sqrt(42),-5/sqrt(42)),...
    complex(-1/sqrt(42),-7/sqrt(42)), complex(-3/sqrt(42),7/sqrt(42)), complex(-3/sqrt(42),5/sqrt(42)),...
    complex(-3/sqrt(42),3/sqrt(42)), complex(-3/sqrt(42),1/sqrt(42)), complex(-3/sqrt(42),-1/sqrt(42)),...
    complex(-3/sqrt(42),-3/sqrt(42)), complex(-3/sqrt(42),-5/sqrt(42)), complex(-3/sqrt(42),-7/sqrt(42)),...
    complex(-5/sqrt(42),7/sqrt(42)), complex(-5/sqrt(42),5/sqrt(42)), complex(-5/sqrt(42),3/sqrt(42)),...
    complex(-5/sqrt(42),1/sqrt(42)), complex(-5/sqrt(42),-1/sqrt(42)), complex(-5/sqrt(42),-3/sqrt(42)),...
    complex(-5/sqrt(42),-5/sqrt(42)), complex(-5/sqrt(42),-7/sqrt(42)), complex(-7/sqrt(42),7/sqrt(42)),...
    complex(-7/sqrt(42),5/sqrt(42)), complex(-7/sqrt(42),3/sqrt(42)), complex(-7/sqrt(42),1/sqrt(42)),...
    complex(-7/sqrt(42),-1/sqrt(42)), complex(-7/sqrt(42),-3/sqrt(42)), complex(-7/sqrt(42),-5/sqrt(42)),...
    complex(-7/sqrt(42),-7/sqrt(42))];

qam_64_comp = zeros(50000, 1);
for i = 1:50000
    sum = 0;
    for j = 1:6
        sum = sum + ranBN((i-1)*6+(j-1)+1);
        sum = sum*2;
    end
    sum = sum/2;
    qam_64(i) = sum;
    qam_64_comp(i) = const_point_qam_64(qam_64(i)+1);
end
%for i = 1:20
    %disp(qam_64_comp(i))
%end

% task4 - transmit over the air

% convert from dBm to watts
n0 = (((10^(-9))*(1e-3))^(1/2));
% disp(n0)

awgn1 = zeros(300000, 1);
for i = 1:300000
    %ori_noise_x = randn();
    %ori_noise_y = randn();
    scaled_noise = (n0) * ori_awgn(i);
    awgn1(i) = scaled_noise;
end

% disp(num2str(awgn(10)))

% disp(num2str(h))

% BPSK
bpsk_yi = zeros(300000, 1);

for i = 1:300000
    bpsk_yi(i) = scaling_h*bpsk_comp(i) + awgn1(i);
end

% QPSK
qpsk_yi = zeros(150000, 1);

for i = 1:150000
    qpsk_yi(i) = scaling_h*qpsk_comp(i) + awgn1(i);
end

% 16-QAM
qam_16_yi = zeros(75000, 1);
for i = 1:75000
    qam_16_yi(i) = scaling_h*qam_16_comp(i) + awgn1(i);
end

% 64-QAM
qam_64_yi = zeros(50000, 1);
for i = 1:50000
    qam_64_yi(i) = scaling_h*qam_64_comp(i) + awgn1(i);
end

%for i = 1:10
    %disp(num2str(qam_64_yi(i)))
%end

% task5 - decode and demodulation

%BPSK
bpsk_xi_prime = zeros(300000, 1);
bpsk_xi_const = zeros(300000, 1);
bpsk_de_rx = zeros(300000, 1);
for i = 1:300000
    bpsk_xi_prime(i) = bpsk_yi(i)/scaling_h;
    tmp_real = real(bpsk_xi_prime(i));
    if tmp_real >= 0
        bpsk_xi_const(i) = const_point_bpsk(1);
        bpsk_de_rx(i) = 0;
    else 
        bpsk_xi_const(i) = const_point_bpsk(2);
        bpsk_de_rx(i) = 1;
    end
end

%for i = 1:10
    %disp(num2str(bpsk_xi_prime(i)))
    %disp(num2str(bpsk_xi_const(i)))
%end

%QPSK
qpsk_xi_prime = zeros(150000, 1);
qpsk_xi_const = zeros(150000, 1);
qpsk_de_rx = zeros(300000, 1);

for i = 1:150000
    qpsk_xi_prime(i) = qpsk_yi(i)/scaling_h;
    tmp_real = real(qpsk_xi_prime(i));
    tmp_imag = imag(qpsk_xi_prime(i));
    if tmp_real >= 0
        if tmp_imag >=0
            qpsk_xi_const(i) = const_point_qpsk(1);
            qpsk_de_rx((i-1)*2+1) = 0;
            qpsk_de_rx((i-1)*2+2) = 0;
        else
            qpsk_xi_const(i) = const_point_qpsk(2);
            qpsk_de_rx((i-1)*2+1) = 0;
            qpsk_de_rx((i-1)*2+2) = 1;
        end
    else
        if tmp_imag >=0
            qpsk_xi_const(i) = const_point_qpsk(3);
            qpsk_de_rx((i-1)*2+1) = 1;
            qpsk_de_rx((i-1)*2+2) = 0;
        else
            qpsk_xi_const(i) = const_point_qpsk(4);
            qpsk_de_rx((i-1)*2+1) = 1;
            qpsk_de_rx((i-1)*2+2) = 1;
        end
    end
end

%for i = 1:5
    %disp(num2str(qpsk_xi_prime(i)))
    %disp(num2str(qpsk_xi_const(i)))
%end

%16-QAM
qam_16_xi_prime = zeros(75000, 1);
qam_16_xi_const = zeros(75000, 1);
qam_16_de_rx = zeros(300000, 1);

for i = 1:75000
    qam_16_xi_prime(i) = qam_16_yi(i)/scaling_h;
    tmp_real1 = real(qam_16_xi_prime(i));
    tmp_imag1 = imag(qam_16_xi_prime(i));
    min_index = 1;
    min_distance = 100;
    for j = 1:16
        tmp_const = const_point_qam_16(j);
        tmp_const_real = real(tmp_const);
        tmp_const_imag = imag(tmp_const);
        tmp_distance = abs(tmp_real1-tmp_const_real) + abs(tmp_imag1-tmp_const_imag);
        if tmp_distance<min_distance
            min_distance = tmp_distance;
            min_index = j;
        end
    end
    qam_16_xi_const(i) = const_point_qam_16(min_index);
    min_index = min_index -1;
    for j = 4:-1:1
        if mod(min_index, 2) == 0
            min_index = min_index/2;
            qam_16_de_rx((i-1)*4+j)=0;
        else
            min_index = (min_index-1)/2;
            qam_16_de_rx((i-1)*4+j)=1;
        end
    end
end

%for i = 10:20
    %disp(num2str(qam_16_xi_const(i)))
%end

%64-QAM
qam_64_xi_prime = zeros(50000, 1);
qam_64_xi_const = zeros(50000, 1);
qam_64_de_rx = zeros(300000, 1);

for i = 1:50000
    qam_64_xi_prime(i) = qam_64_yi(i)/scaling_h;
    tmp_real1 = real(qam_64_xi_prime(i));
    tmp_imag1 = imag(qam_64_xi_prime(i));
    min_index = 1;
    min_distance = 100;
    for j = 1:64
        tmp_const = const_point_qam_64(j);
        tmp_const_real = real(tmp_const);
        tmp_const_imag = imag(tmp_const);
        tmp_distance = abs(tmp_real1-tmp_const_real) + abs(tmp_imag1-tmp_const_imag);
        if tmp_distance<min_distance
            min_distance = tmp_distance;
            min_index = j;
        end
    end
    qam_64_xi_const(i) = const_point_qam_64(min_index);
    min_index = min_index - 1;
    for j = 6:-1:1
        if mod(min_index, 2)==0
            min_index = min_index/2;
            qam_64_de_rx((i-1)*6+j) = 0;
        else
            min_index = (min_index-1)/2;
            qam_64_de_rx((i-1)*6+j) = 1;
        end
    end
end

%for i = 10:20
    %disp(num2str(qam_64_xi_const(i)))
%end

%task6 - calculate SNR

%BPSK
bpsk_emp_error = zeros(300000, 1);
total_error = 0;

qpsk_emp_error = zeros(150000, 1);
qam_16_emp_error = zeros(75000, 1);
qam_64_emp_error = zeros(50000, 1);

for i = 1:300000
    bpsk_emp_error(i) = bpsk_xi_prime(i) - bpsk_xi_const(i);
    total_error = total_error + (real(bpsk_emp_error(i)))^2 + (imag(bpsk_emp_error(i)))^2;
end
for i = 1:150000
    qpsk_emp_error(i) = qpsk_xi_prime(i) - qpsk_xi_const(i);
    total_error = total_error + (real(qpsk_emp_error(i)))^2 + (imag(qpsk_emp_error(i)))^2;
end
for i = 1:75000
    qam_16_emp_error(i) = qam_16_xi_prime(i) - qam_16_xi_const(i);
    total_error = total_error + (real(qam_16_emp_error(i)))^2 + (imag(qam_16_emp_error(i)))^2;
end
for i = 1:50000
    qam_64_emp_error(i) = qam_64_xi_prime(i) - qam_64_xi_const(i);
    total_error = total_error + (real(qam_64_emp_error(i)))^2 + (imag(qam_64_emp_error(i)))^2;
end

bpsk_avg_error_watt = total_error/(300000+150000+75000+50000);
scaling_h_power = (real(scaling_h))^2 + (imag(scaling_h))^2;
bpsk_avg_noise_emp = (bpsk_avg_error_watt*scaling_h_power);
t6_noise_emp_watt = ['Empirical average noise power in Watt: ', num2str(bpsk_avg_noise_emp)];
disp(t6_noise_emp_watt)

thr_noise_watt = (10^(-90/10))*(1e-3);
t6_noise_the_watt = ['Theoretical average noise power in Watt: ', num2str(thr_noise_watt)];
disp(t6_noise_the_watt)
%disp(bpsk_avg_noise_emp)
% disp(scaling_h_power)

bpsk_avg_noise_dbm = 10*log10(bpsk_avg_noise_emp/1e-3);
t6_noise_emp_dbm = ['Empirical average noise power in dBm: ', num2str(bpsk_avg_noise_dbm)];
disp(t6_noise_emp_dbm)
%disp(bpsk_avg_noise_dbm)
% bpsk_snr_emp = rx_power_watt/bpsk_avg_noise_emp;
thr_noise_dbm = -90;
t6_noise_the_dbm = ['Theoretical average noise power in dBm: ', num2str(thr_noise_dbm)];
disp(t6_noise_the_dbm)

bpsk_snr_db_emp = 10*log10(rx_power_watt/1e-3) - bpsk_avg_noise_dbm;
t6_bpsk1 = ['Empirical SNRdB: ', num2str(bpsk_snr_db_emp)];
disp(t6_bpsk1)

bpsk_snr_db_the = rx_power_dbm - (-90);
t6_bpsk = ['Theoretical SNRdB: ', num2str(bpsk_snr_db_the)];
disp(t6_bpsk)

%task7 - calculate throughput

%BPSK
%BER
bpsk_total_bit_error = 0;
for i = 1:300000
    if bpsk_de_rx(i) ~= bpsk(i)
        bpsk_total_bit_error = bpsk_total_bit_error +1;
    end
end

for i = 1:10
    %disp(bpsk_de_rx(i))
    %disp(bpsk(i))
end
bpsk_ber = bpsk_total_bit_error/300000;
%disp(bpsk_total_bit_error)
t7_bpsk_ber = ['Empirical BER of BPSK: ', num2str(bpsk_ber)];
disp(t7_bpsk_ber)
%PDR
bpsk_pdr_the = (1-bpsk_ber)^(500*8);
%disp(bpsk_pdr_the)

bpsk_packet_total_success = 0;
for i = 1:75
    bpsk_sign = 1;
    for j = 1:4000
        if bpsk_de_rx((i-1)*4000+(j-1)+1) ~= bpsk((i-1)*4000+(j-1)+1)
            bpsk_sign = 0;
        end
    end
    if bpsk_sign == 1
        bpsk_packet_total_success = bpsk_packet_total_success + 1;
    end
end
%disp(bpsk_packet_total_success)
bpsk_pdr_emp = bpsk_packet_total_success/75;
%disp(bpsk_pdr_emp)
%throughput
bpsk_throughput_emp = 0;
bpsk_throughput_emp = bpsk_pdr_emp*((1/3.2)*1e6);
y1(count) = bpsk_throughput_emp;
%disp(y1(count))
%disp(count)
%disp(bpsk_throughput_emp)
t7_bpsk_throughput = ['Empirical Throughput of BPSK(bps): ', num2str(bpsk_throughput_emp)];
disp(t7_bpsk_throughput)

bpsk_throughput_thr = 0;
bpsk_throughput_thr = bpsk_pdr_the*((1/3.2)*1e6);
%disp(bpsk_throughput_thr)
t7_bpsk1_throughput = ['Theoretical Throughput of BPSK(bps): ', num2str(bpsk_throughput_thr)];
disp(t7_bpsk1_throughput)

%QPSK
%BER
qpsk_total_bit_error = 0;
for i = 1:300000
    if qpsk_de_rx(i) ~= bpsk(i)
        qpsk_total_bit_error = qpsk_total_bit_error +1;
    end
end

qpsk_ber = qpsk_total_bit_error/300000;
t7_qpsk_ber = ['Empirical BER of QPSK: ', num2str(qpsk_ber)];
disp(t7_qpsk_ber)
%PDR
qpsk_pdr_the = (1-qpsk_ber)^(500*8);
%disp(qpsk_pdr_the)

qpsk_packet_total_success = 0;
for i = 1:75
    qpsk_sign = 1;
    for j = 1:4000
        if qpsk_de_rx((i-1)*4000+(j-1)+1) ~= bpsk((i-1)*4000+(j-1)+1)
            qpsk_sign = 0;
        end
    end
    if qpsk_sign == 1
        qpsk_packet_total_success = qpsk_packet_total_success + 1;
    end
end
%disp(qpsk_packet_total_success)
qpsk_pdr_emp = qpsk_packet_total_success/75;
%disp(qpsk_pdr_emp)

%throughput
qpsk_throughput_emp = 0;
qpsk_throughput_emp = qpsk_pdr_emp*((2/3.2)*1e6);
y2(distance/50) = qpsk_throughput_emp;
t7_qpsk_throughput = ['Empirical Throughput of QPSK(bps): ', num2str(qpsk_throughput_emp)];
disp(t7_qpsk_throughput)

qpsk_throughput_thr = 0;
qpsk_throughput_thr = qpsk_pdr_the*((2/3.2)*1e6);
%disp(qpsk_throughput_thr)
t7_qpsk1_throughput = ['Theoretical Throughput of QPSK(bps): ', num2str(qpsk_throughput_thr)];
disp(t7_qpsk1_throughput)

%16-QAM
qam_16_total_bit_error = 0;
for i = 1:300000
    if qam_16_de_rx(i) ~= bpsk(i)
        qam_16_total_bit_error = qam_16_total_bit_error +1;
    end
end

for i = 1:10
    %disp(qam_16_de_rx(i))
    %disp(bpsk(i))
end
qam_16_ber = qam_16_total_bit_error/300000;
t7_qam16_ber = ['Empirical BER of 16-QAM: ', num2str(qam_16_ber)];
disp(t7_qam16_ber)
%PDR
qam_16_pdr_the = (1-qam_16_ber)^(500*8);
%disp(qam_16_pdr_the)

qam_16_packet_total_success = 0;
for i = 1:75
    qam_16_sign = 1;
    for j = 1:4000
        if qam_16_de_rx((i-1)*4000+(j-1)+1) ~= bpsk((i-1)*4000+(j-1)+1)
            qam_16_sign = 0;
        end
    end
    if qam_16_sign == 1
        qam_16_packet_total_success = qam_16_packet_total_success + 1;
    end
end
%disp(qpsk_packet_total_success)
qam_16_pdr_emp = qam_16_packet_total_success/75;
%disp(qpsk_pdr_emp)

%throughput
qam_16_throughput_emp = 0;
qam_16_throughput_emp = qam_16_pdr_emp*((4/3.2)*1e6);
y3(distance/50) = qam_16_throughput_emp;
t7_qam_16_throughput = ['Empirical Throughput of 16-QAM(bps): ', num2str(qam_16_throughput_emp)];
disp(t7_qam_16_throughput)

qam_16_throughput_thr = 0;
qam_16_throughput_thr = qam_16_pdr_the*((4/3.2)*1e6);
%disp(qam_16_throughput_thr)
t7_qam1_16_throughput = ['Theoretical Throughput of 16-QAM(bps): ', num2str(qam_16_throughput_thr)];
disp(t7_qam1_16_throughput)

%64-QAM
qam_64_total_bit_error = 0;
for i = 1:300000
    if qam_64_de_rx(i) ~= bpsk(i)
        qam_64_total_bit_error = qam_64_total_bit_error +1;
    end
end

for i = 1:10
    %disp(qam_64_de_rx(i))
    %disp(bpsk(i))
end
qam_64_ber = qam_64_total_bit_error/300000;
t7_qam64_ber = ['Empirical BER of 64-QAM: ', num2str(qam_64_ber)];
disp(t7_qam64_ber)
%PDR
qam_64_pdr_the = (1-qam_64_ber)^(500*8);
%disp(qam_64_pdr_the)

qam_64_packet_total_success = 0;
for i = 1:75
    qam_64_sign = 1;
    for j = 1:4000
        if qam_64_de_rx((i-1)*4000+(j-1)+1) ~= bpsk((i-1)*4000+(j-1)+1)
            qam_64_sign = 0;
        end
    end
    if qam_64_sign == 1
        qam_64_packet_total_success = qam_64_packet_total_success + 1;
    end
end
%disp(qpsk_packet_total_success)
qam_64_pdr_emp = qam_64_packet_total_success/75;
%disp(qpsk_pdr_emp)

%throughput
qam_64_throughput_emp = 0;
qam_64_throughput_emp = qam_64_pdr_emp*((6/3.2)*1e6);
y4(distance/50) = qam_64_throughput_emp;
t7_qam_64_throughput = ['Empirical Throughput of 64-QAM(bps): ', num2str(qam_64_throughput_emp)];
disp(t7_qam_64_throughput)

qam_64_throughput_thr = 0;
qam_64_throughput_thr = qam_64_pdr_the*((6/3.2)*1e6);
%disp(qam_16_throughput_thr)
t7_qam1_64_throughput = ['Theoretical Throughput of 64-QAM(bps): ', num2str(qam_64_throughput_thr)];
disp(t7_qam1_64_throughput)

%optimal modulation scheme
max_throughput = bpsk_throughput_emp;
max_sel = 1;

if qpsk_throughput_emp>max_throughput
    max_throughput = qpsk_throughput_emp;
    max_sel = 2;
end
if qam_16_throughput_emp>max_throughput
    max_throughput = qam_16_throughput_emp;
    max_sel = 3;
end
if qam_64_throughput_emp>max_throughput
    max_throughput = qam_64_throughput_emp;
    max_sel = 4;
end

if bpsk_throughput_emp == 0
    t7_optimal_sel = ['Optimal modulation scheme: None'];
elseif max_sel == 1
    t7_optimal_sel = ['Optimal modulation scheme: BPSK'];
elseif max_sel == 2
    t7_optimal_sel = ['Optimal modulation scheme: QPSK'];
elseif max_sel == 3
    t7_optimal_sel = ['Optimal modulation scheme: 16-QAM'];
else
    t7_optimal_sel = ['Optimal modulation scheme: 64-QAM'];
end

disp(t7_optimal_sel)

%task8 - plot figures

%Throughput v.s. distance
count = count + 1;
end

dist = 50:50:600;

plot(dist, y1, 'g', dist, y2, 'r', dist, y3, 'b', dist, y4, 'k')
xlabel('distance')
ylabel('throughput')
legend('BPSK', 'QPSK', '16-QAM', '64-QAM')

tx_gain = 1;
rx_gain = 1;
tx_power_dbm = 10;
wavelength = 3e+8/2.4e+9;
%constellation diagram
%BPSK
bpsk_dist = 40;
rx_power_dbm = tx_power_dbm + 10*log10(tx_gain) + 10*log10(rx_gain) + 20*log10(wavelength/(4*pi*bpsk_dist));
rx_power_watt = (10^(rx_power_dbm/10))*(1e-3);
h = complex(a, b);
scaling_h = h * ((rx_power_watt)^(1/2));
n0 = (((10^(-9))*(1e-3))^(1/2));
const_point_bpsk = [complex(1, 0), complex(-1, 0)];
bpsk_comp = zeros(300000, 1);
awgn1 = zeros(300000, 1);
bpsk = zeros(300000, 1);
bpsk_yi = zeros(300000, 1);

for i = 1:300000
    bpsk(i) = ranBN(i);
    bpsk_comp(i) = const_point_bpsk(bpsk(i)+1);
end
for i = 1:300000
    %ori_noise_x = randn();
    %ori_noise_y = randn();
    scaled_noise = (n0) * ori_awgn(i);
    awgn1(i) = scaled_noise;
end
for i = 1:300000
    bpsk_yi(i) = scaling_h*bpsk_comp(i) + awgn1(i);
end
for i = 1:300000
    bpsk_xi_prime(i) = bpsk_yi(i)/scaling_h;
end

bpsk_real_part = real(bpsk_xi_prime);
bpsk_imag_part = imag(bpsk_xi_prime);

figure;
scatter(bpsk_real_part, bpsk_imag_part, 'filled', 'SizeData',1);
xlabel('I');
ylabel('Q');
title('BPSK constellation diagram');
grid on;

%QPSK
qpsk_dist = 20;
rx_power_dbm = tx_power_dbm + 10*log10(tx_gain) + 10*log10(rx_gain) + 20*log10(wavelength/(4*pi*qpsk_dist));
rx_power_watt = (10^(rx_power_dbm/10))*(1e-3);
h = complex(a, b);
scaling_h = h * ((rx_power_watt)^(1/2));
n0 = (((10^(-9))*(1e-3))^(1/2));
const_point_qpsk = [complex(1/sqrt(2), 1/sqrt(2)), complex(1/sqrt(2), -1/sqrt(2)),...
complex(-1/sqrt(2), 1/sqrt(2)), complex(-1/sqrt(2), -1/sqrt(2))];

qpsk_comp = zeros(150000, 1);
qpsk = zeros(150000, 1);
qpsk_yi = zeros(150000, 1);

for i = 1:150000
    sum = 0;
    for j = 1:2
        sum = sum + ranBN((i-1)*2+(j-1)+1);
        sum = sum * 2;
    end
    sum = sum/2;
    qpsk(i) = sum;
    qpsk_comp(i) = const_point_qpsk(qpsk(i)+1);
    % disp(num2str(ranBN(i)))
end
for i = 1:150000
    qpsk_yi(i) = scaling_h*qpsk_comp(i) + awgn1(i);
end
for i = 1:150000
    qpsk_xi_prime(i) = qpsk_yi(i)/scaling_h;
end

qpsk_real_part = real(qpsk_xi_prime);
qpsk_imag_part = imag(qpsk_xi_prime);

figure;
scatter(qpsk_real_part, qpsk_imag_part, 'filled', 'SizeData',1);
xlabel('I');
ylabel('Q');
title('QPSK constellation diagram');
grid on;

%16-QAM
qam_16_dist = 10;
rx_power_dbm = tx_power_dbm + 10*log10(tx_gain) + 10*log10(rx_gain) + 20*log10(wavelength/(4*pi*qam_16_dist));
rx_power_watt = (10^(rx_power_dbm/10))*(1e-3);
h = complex(a, b);
scaling_h = h * ((rx_power_watt)^(1/2));
n0 = (((10^(-9))*(1e-3))^(1/2));
const_point_qam_16 = [complex(3/sqrt(10),3/sqrt(10)), complex(3/sqrt(10),1/sqrt(10)), complex(3/sqrt(10), -1/sqrt(10)),... 
    complex(3/sqrt(10),-3/sqrt(10)), complex(1/sqrt(10),3/sqrt(10)), complex(1/sqrt(10),1/sqrt(10)),...
    complex(1/sqrt(10),-1/sqrt(10)), complex(1/sqrt(10),-3/sqrt(10)), complex(-1/sqrt(10),3/sqrt(10)),...
    complex(-1/sqrt(10),1/sqrt(10)), complex(-1/sqrt(10),-1/sqrt(10)), complex(-1/sqrt(10),-3/sqrt(10)),...
    complex(-3/sqrt(10),3/sqrt(10)), complex(-3/sqrt(10),1/sqrt(10)), complex(-3/sqrt(10),-1/sqrt(10)),...
    complex(-3/sqrt(10),-3/sqrt(10))];

qam_16 = zeros(75000, 1);
qam_16_yi = zeros(75000, 1);
qam_16_comp = zeros(75000, 1);

for i = 1:75000
    sum = 0;
    for j = 1:4
        sum = sum + ranBN((i-1)*4+(j-1)+1);
        sum = sum*2;
    end
    sum = sum/2;
    qam_16(i) = sum;
    qam_16_comp(i) = const_point_qam_16(qam_16(i)+1);
end
for i = 1:75000
    qam_16_yi(i) = scaling_h*qam_16_comp(i) + awgn1(i);
end
for i = 1:75000
    qam_16_xi_prime(i) = qam_16_yi(i)/scaling_h;
end
qam_16_real_part = real(qam_16_xi_prime);
qam_16_imag_part = imag(qam_16_xi_prime);

figure;
scatter(qam_16_real_part, qam_16_imag_part, 'filled', 'SizeData',1);
xlabel('I');
ylabel('Q');
title('16-QAM constellation diagram');
grid on;

%64-QAM
qam_64_dist = 5;
rx_power_dbm = tx_power_dbm + 10*log10(tx_gain) + 10*log10(rx_gain) + 20*log10(wavelength/(4*pi*qam_64_dist));
rx_power_watt = (10^(rx_power_dbm/10))*(1e-3);
h = complex(a, b);
scaling_h = h * ((rx_power_watt)^(1/2));
n0 = (((10^(-9))*(1e-3))^(1/2));
const_point_qam_64 = [complex(7/sqrt(42),7/sqrt(42)), complex(7/sqrt(42),5/sqrt(42)), complex(7/sqrt(42),3/sqrt(42)),...
    complex(7/sqrt(42),1/sqrt(42)), complex(7/sqrt(42),-1/sqrt(42)), complex(7/sqrt(42),-3/sqrt(42)),...
    complex(7/sqrt(42),-5/sqrt(42)), complex(7/sqrt(42),-7/sqrt(42)), complex(5/sqrt(42),7/sqrt(42)),...
    complex(5/sqrt(42),5/sqrt(42)), complex(5/sqrt(42),3/sqrt(42)), complex(5/sqrt(42),1/sqrt(42)),...
    complex(5/sqrt(42),-1/sqrt(42)), complex(5/sqrt(42),-3/sqrt(42)), complex(5/sqrt(42),-5/sqrt(42)),...
    complex(5/sqrt(42),-7/sqrt(42)), complex(3/sqrt(42),7/sqrt(42)), complex(3/sqrt(42),5/sqrt(42)),...
    complex(3/sqrt(42),3/sqrt(42)), complex(3/sqrt(42),1/sqrt(42)), complex(3/sqrt(42),-1/sqrt(42)),...
    complex(3/sqrt(42),-3/sqrt(42)), complex(3/sqrt(42),-5/sqrt(42)), complex(3/sqrt(42),-7/sqrt(42)),...
    complex(1/sqrt(42),7/sqrt(42)), complex(1/sqrt(42),5/sqrt(42)), complex(1/sqrt(42),3/sqrt(42)),...
    complex(1/sqrt(42),1/sqrt(42)), complex(1/sqrt(42),-1/sqrt(42)), complex(1/sqrt(42),-3/sqrt(42)),...
    complex(1/sqrt(42),-5/sqrt(42)), complex(1/sqrt(42),-7/sqrt(42)), complex(-1/sqrt(42),7/sqrt(42)),...
    complex(-1/sqrt(42),5/sqrt(42)), complex(-1/sqrt(42),3/sqrt(42)), complex(-1/sqrt(42),1/sqrt(42)),...
    complex(-1/sqrt(42),-1/sqrt(42)), complex(-1/sqrt(42),-3/sqrt(42)), complex(-1/sqrt(42),-5/sqrt(42)),...
    complex(-1/sqrt(42),-7/sqrt(42)), complex(-3/sqrt(42),7/sqrt(42)), complex(-3/sqrt(42),5/sqrt(42)),...
    complex(-3/sqrt(42),3/sqrt(42)), complex(-3/sqrt(42),1/sqrt(42)), complex(-3/sqrt(42),-1/sqrt(42)),...
    complex(-3/sqrt(42),-3/sqrt(42)), complex(-3/sqrt(42),-5/sqrt(42)), complex(-3/sqrt(42),-7/sqrt(42)),...
    complex(-5/sqrt(42),7/sqrt(42)), complex(-5/sqrt(42),5/sqrt(42)), complex(-5/sqrt(42),3/sqrt(42)),...
    complex(-5/sqrt(42),1/sqrt(42)), complex(-5/sqrt(42),-1/sqrt(42)), complex(-5/sqrt(42),-3/sqrt(42)),...
    complex(-5/sqrt(42),-5/sqrt(42)), complex(-5/sqrt(42),-7/sqrt(42)), complex(-7/sqrt(42),7/sqrt(42)),...
    complex(-7/sqrt(42),5/sqrt(42)), complex(-7/sqrt(42),3/sqrt(42)), complex(-7/sqrt(42),1/sqrt(42)),...
    complex(-7/sqrt(42),-1/sqrt(42)), complex(-7/sqrt(42),-3/sqrt(42)), complex(-7/sqrt(42),-5/sqrt(42)),...
    complex(-7/sqrt(42),-7/sqrt(42))];

qam_64 = zeros(50000, 1);
qam_64_yi = zeros(50000, 1);
qam_64_comp = zeros(50000, 1);

for i = 1:50000
    sum = 0;
    for j = 1:6
        sum = sum + ranBN((i-1)*6+(j-1)+1);
        sum = sum*2;
    end
    sum = sum/2;
    qam_64(i) = sum;
    qam_64_comp(i) = const_point_qam_64(qam_64(i)+1);
end
for i = 1:50000
    qam_64_yi(i) = scaling_h*qam_64_comp(i) + awgn1(i);
end
for i = 1:50000
    qam_64_xi_prime(i) = qam_64_yi(i)/scaling_h;
end
qam_64_real_part = real(qam_64_xi_prime);
qam_64_imag_part = imag(qam_64_xi_prime);

figure;
scatter(qam_64_real_part, qam_64_imag_part, 'filled', 'SizeData',1);
xlabel('I');
ylabel('Q');
title('64-QAM constellation diagram');
grid on;