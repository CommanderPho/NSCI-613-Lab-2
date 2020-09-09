% Produce Final Results: Lab 2
% After calling run_HodHux:

spikeFrequency = spikeFrequency_mean;
[min_spikeFrequency, min_freq_i] = min(spikeFrequency);
[max_spikeFrequency, max_freq_i] = max(spikeFrequency);
min_current_x = pulsei_var(min_freq_i);
max_current_x = pulsei_var(max_freq_i);
fprintf('Min: (%.6g, %.6g)\n',min_current_x,min_spikeFrequency);
fprintf('Max: (%.6g, %.6g)\n',max_current_x,max_spikeFrequency);
% [1.1, 49.49]
% 
% Compute the f_{c}:
f_c = min_spikeFrequency; 


figure(2)
clf
plot(pulsei_var',spikeFrequency)
if is_modified_potassium_curve_mode
	title('f-I relation (Modified K+ Mode) for I_{app} between 0 and 60 [\mu A / cm^{2}]')
else
	title('f-I relation for I_{app} between 0 and 60 [\mu A / cm^{2}]')
end
xlabel('Applied Depolarizing Current [\mu A / cm^{2}]')
ylabel('Firing rate [Hz]')
xlim([0,60])
ylim([0,max_spikeFrequency*1.1])
hold on;
scatter([min_current_x, max_current_x], [min_spikeFrequency, max_spikeFrequency],'r')


%% Multi-voltage curve plots
% Multi-subplot version:
figure(4)
clf
hold off;

if is_modified_potassium_curve_mode
	searchPlotValues = [1.0, 1.1, 1.2];
else
% 	searchPlotValues = [0.0, 0.1, 1.0];
% 	searchPlotValues = [2.1, 2.2, 2.3, 2.4];
% 	searchPlotValues = [6.0, 6.1, 6.2, 6.3, 6.4, 6.5];
% 	searchPlotValues = [61.0, 70.0, 82.0, 83.0, 170, 172, 173];

end


interestingPlotIndicies = zeros(length(searchPlotValues),1);
for i=1:length(searchPlotValues)
	found_indicies = find(abs(resultsTable.AppliedCurrent-searchPlotValues(i)) < 0.001);
	interestingPlotIndicies(i) = found_indicies;
% 	interestingPlotIndicies(i) = find(resultsTable.AppliedCurrent == searchPlotValues(i));
end

% interestingPlotIndicies = find(resultsTable.AppliedCurrent == searchPlotValues);

%interestingPlotIndicies = [41, 42, 43, 44, 45]; % P2a: 6mA, 6.1mA, 6.2mA, 6.3mA, 6.4mA
% interestingPlotIndicies = [41]; % P2b: 6mA, Subthreshold oscillations.
% interestingPlotIndicies = [44]; % P2b: 6.3mA, Subthreshold oscillations.


%interestingPlotIndicies = [2.2, 2.3];
%interestingPlotIndicies = [3, 4, 5]; % Single-spike firing transition

%interestingPlotIndicies = [3]; % For looking a subthreshold oscillations.

%interestingPlotIndicies = [41, 42, 43, 44, 45, 46]; % Continuous Transition

%interestingPlotIndicies = [53, 73]; % P3: 60mA, 80mA 
% interestingPlotIndicies = [74]; % P3: 81mA, 82mA 
%interestingPlotIndicies = [13, 33, 53, 73, 74]; % P3: 20mA, 40mA, 60mA, 80mA, 81mA
%interestingPlotIndicies = [3, 4, 5, 13, 33, 53, 73, 74]; % P3: 20mA, 40mA, 60mA, 80mA, 81mA
% interestingPlotIndicies = [101]; % P3: 108mA

% interestingPlotIndicies = [581, 781, 791]; % P3: 60mA, 80mA, 81mA

% interestingPlotIndicies = [10, 11, 12]; % P4: 1mA, 1.1mA, ...

% interestingPlotIndicies = [172, 173, 180]; % P4: 106mA, 107mA, 108mA ...

num_active_indices = length(interestingPlotIndicies);
curr_voltage_strings = {};
curr_sub_strings = {};

curr_subthresholdFrequency_mean = zeros(num_active_indices,1);
curr_subthresholdFrequency_last = zeros(num_active_indices,1);

% Loop through the results
for i=1:num_active_indices
	active_i = interestingPlotIndicies(i);
	active_applied_current_value = resultsTable.AppliedCurrent(active_i);
	active_applied_spikeCount_value = resultsTable.spikeCounts(active_i);
	active_applied_spikeFrequency_value = resultsTable.spikeFrequency_mean(active_i);
	curr_time_t_data = time_t{active_i};
	curr_Vm_data = voltageTraces{active_i};
	
	%% Test Findpeaks:
	% determine spike times and interspike intervals
	[curr_all_peaks, curr_all_peaks_indxs] = findpeaks(curr_Vm_data,'MinPeakProminence',0.5); % Get the subthreshold peaks only
	
	curr_peak_is_subthreshold = (curr_all_peaks < -10);
	curr_subthreshold_peak_indices = curr_all_peaks_indxs(curr_peak_is_subthreshold);
	curr_subthreshold_peaks = curr_all_peaks(curr_peak_is_subthreshold);

% 	% All Peaks mode:
% 	curr_indxs = curr_all_peaks_indxs;
% 	curr_peaks = curr_all_peaks;

	
	% Subthreshold Only mode:
	curr_indxs = curr_subthreshold_peak_indices;
	curr_peaks = curr_subthreshold_peaks;
	
% 	curr_peakCounts = length(curr_indxs);
	
	curr_subthreshold_times = curr_time_t_data(curr_indxs);
	curr_subthreshold_intervals = diff(curr_subthreshold_times);
	
	if ~isempty(curr_subthreshold_intervals)
		last_IPI_seconds = curr_subthreshold_intervals(end) / 1000;
		curr_subthresholdFrequency_last(i) = 1 ./ last_IPI_seconds;
		
		mean_ISI_seconds = mean(curr_subthreshold_intervals) / 1000;
		curr_subthresholdFrequency_mean(i) = 1 ./ mean_ISI_seconds;
	else
		curr_subthresholdFrequency_last(i) = NaN;
		curr_subthresholdFrequency_mean(i) = NaN;
	end
	
	
	% Plot the current data
	curr_subplot_handle = subplot(num_active_indices,1,i);
% 	curr_plot_handle = plot(curr_time_t_data, curr_Vm_data);
	curr_plot_handle = plot(curr_time_t_data, curr_Vm_data, curr_subthreshold_times, curr_peaks, '*');
	
	
	%hold on;
	curr_voltage_string = sprintf('%f [\\mu A / cm^{2}]', active_applied_current_value);
	curr_sub_string = sprintf('Spikes: %d, SpikeFreq: %.6g [Hz]', active_applied_spikeCount_value, active_applied_spikeFrequency_value);
	curr_sub_threshold_string = sprintf('Subthreshold Freq: %.6g [Hz]', curr_subthresholdFrequency_mean(i));
	
	curr_sub_string = [curr_sub_string, ' ', curr_sub_threshold_string];
	
	curr_voltage_strings = [curr_voltage_strings; curr_voltage_string];
	curr_sub_strings = [curr_sub_strings; curr_sub_string];
	
	ylabel('Voltage [mV]')
	title(curr_voltage_string)
	fprintf('curr_voltage_string[%d]: %s\n',active_i, curr_voltage_string);
	fprintf('\t curr_sub_string[%d]: %s\n',active_i, curr_sub_string);
	
end
xlabel('Time [ms]')
% ylabel('Voltage [mV]')
% legend(curr_voltage_strings)
%title('Depolariziation Block Voltage Plot')
if is_modified_potassium_curve_mode
	sgtitle('I_{app} threshold determination for (Modified K+ Mode)')
else
	sgtitle('I_{app} threshold determination')
end




