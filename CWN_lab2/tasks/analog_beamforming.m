function [pathloss_rx_dbm, rx1_SNR_dbm, rx2_INR_dbm] = analog_beamforming(codebook_size, P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, tx_antenna_number)
addpath ./ewa_function;

% Generate codebook for analog beamforming
tx_beam_direction = 0:codebook_size:180; % degree
d = 0.5;                      % Distance between antennas (multiple number of wavelength)
% Generate gain table: calculate power gain of beam at each direction(beam direction * sector angle)
% gain = gain(beam, sector angle)
resolution = 360;                                           % Number of angles dividing 180 degrees
a = zeros(numel(tx_beam_direction), tx_antenna_number);     % Antennas' array coefficient(weights), direction * number of antenna
A = zeros(numel(tx_beam_direction), resolution);            % Antennas' array factor
gain = zeros(numel(tx_beam_direction), resolution);         % Antennas' power gain
phi = (1 : resolution) * pi / resolution;                   % Equally-spaced over [0, pi] into resolution angle
phi_deg = phi * 180 / pi;                                   % Convert phi to degree
psi = 2 * pi * d * cos(phi);                                % Antenna phase shift of angle phi

for tx_beam_idx = 1:numel(tx_beam_direction)
    ph0 = tx_beam_direction(tx_beam_idx);
    a(tx_beam_idx,:) = uniform(d, ph0, tx_antenna_number);  % beam steering weights
    A(tx_beam_idx,:) = dtft(a(tx_beam_idx,:), -psi);        % Array factor, note dtft(a,psi)=dtft(a,-psi)
    gain(tx_beam_idx,:) = abs(A(tx_beam_idx,:)).^2;         % Power gain
end

% TODO1: Calculate the closest_sector_index and beam_index using gain table
% Hint: rx1_sector_index: Find the closest sector angle based on rx1_theta and rx2_theta
% Hint: rx1_beam_index: Find the optimal beam with the highest power gain for rx1
rx1_distance = norm(rx1_location - tx_location); % Distance between receiver 1 and transmitter (m)
rx2_distance = norm(rx2_location - tx_location); % Distance between receiver 2 and transmitter (m)
rx1_theta = atan2(rx1_location(2), rx1_location(1)); % Direction of receiver 1 (rad)
rx2_theta = atan2(rx2_location(2), rx2_location(1)); % Direction of receiver 2 (rad)
%disp(gain)
%disp(rx1_theta)
%disp(rx2_theta)
%disp(rx1_location)
%disp(rx2_location)
% calculate closest index(rx1)
rx1_sector_index = 0;
rx1_smallest_dif = Inf;
for i = 1:1:360
    if abs(phi(i) - rx1_theta) < rx1_smallest_dif
        rx1_smallest_dif = abs(phi(i) - rx1_theta);
        rx1_sector_index = i;
    end
    %disp(rx1_smallest_dif)
    %disp(rx1_sector_index)
end

% rx2_sector_index
rx2_sector_index = 0;
rx2_smallest_dif = Inf;
for i = 1:1:360
    if abs(phi(i) - rx2_theta) < rx2_smallest_dif
        rx2_smallest_dif = abs(phi(i) - rx2_theta);
        rx2_sector_index = i;
    end
    %disp(rx2_smallest_dif)
    %disp(rx2_sector_index)
end
%disp(rx2_sector_index)

% calculate optimal beam index
rx1_beam_index = 0;
rx1_largest_gain = 0;
for tx_beam_idx = 1:numel(tx_beam_direction)
    if gain(tx_beam_idx, rx1_sector_index) > rx1_largest_gain
        rx1_largest_gain = gain(tx_beam_idx, rx1_sector_index);
        rx1_beam_index = tx_beam_idx;
    end
end
%disp(rx1_beam_index)

% Uncomment to plot the beam pattern and show the results of optimal beam
%figure(1); 
%polarplot(phi, gain(rx1_beam_index,:));
%title('Power Gain in different directions');
% 
%figure(2); 
%plot(phi_deg, gain(rx1_beam_index, :)); 
%xlabel('Angle (degree)');
%xlim([0 180]);
%ylabel('Power Gain');
%title('Power Gain vs. Angle');
% 
%fprintf('Beam Scanning Results\n');
%fprintf('\tReceiver 1 Distance: %f m\tAngle of departure: %2f degree\n', rx1_distance, phi_deg(rx1_sector_index));
%fprintf('\tReceiver 2 Distance: %f m\tAngle of departure: %2f degree\n', rx2_distance, phi_deg(rx2_sector_index));
%fprintf('\tClosest beam direction is %.2f (index: %d)\n', tx_beam_direction(rx1_beam_index), rx1_beam_index);
%fprintf('\t%d-th beam is the best beam for rx1 (with power gain %f)\n', rx1_beam_index, gain(rx1_beam_index, rx1_sector_index));
%fprintf('\t%d-th beam side lobe for rx2 (with power gain %f)\n', rx1_beam_index, gain(rx1_beam_index, rx2_sector_index));

% TODO2: Calculate receiving power and SNR (Prx-noise)
% Hint: Use friis equation with gain
freq = 24e9;
rx1_gain = 0;
rx2_gain = 0;
rx1_SNR_dbm = 0;
rx2_INR_dbm = 0;

% calculate Gtx of rx1 and rx2
rx1_gain = gain(rx1_beam_index, rx1_sector_index);
rx2_gain = gain(rx1_beam_index, rx2_sector_index);

% calculate Ptx of rx1 and rx2 by friis-equation
pathloss_rx_dbm = friis_equation(freq, rx1_gain, 1, rx1_distance) + P_tx_dBm;
pathloss_interf_dbm = friis_equation(freq, rx2_gain, 1, rx2_distance) + P_tx_dBm;

%disp(pathloss_rx_dbm)
%disp(pathloss_interf_dbm)
%disp(N0_dBm)

% calculate SNR and INR by SNRdb = Psignal,dbm - Pnoise,dbm
rx1_SNR_dbm = (pathloss_rx_dbm)-(N0_dBm);
rx2_INR_dbm = (pathloss_interf_dbm)-(N0_dBm);

