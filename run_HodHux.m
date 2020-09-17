% script to numerically simulate the HH model
% calls function HHEquations.m which must be in same directory as this file

% set simulation timespan in ms
tspan = [0 1100];


ComputationalOptions.shouldComputePeaks = false;


% Set HH Values:
HH_Parameters.c=1;
HH_Parameters.gNa=120;
HH_Parameters.gK=36;
HH_Parameters.gl=0.3;
HH_Parameters.eNa=50;
HH_Parameters.eK=-77;
HH_Parameters.el=-54.4;

			
% set initial conditions [V(t=0); m(t=0); h(t=0); n(t=0)]
v0 = -64.9997;
m0 = 0.0529;
h0 = 0.5961;
n0 = 0.3176;
ICs = [v0; m0; h0; n0];


pulsei_var = 1:0.1:20;  % REDO: All

num_iterations = length(pulsei_var);

is_modified_potassium_curve_mode = false;
if is_modified_potassium_curve_mode
	curr_mode_string = 'Modified K+ Mode';
else
	curr_mode_string = '';	
end

if is_modified_potassium_curve_mode
	% Delayed Potassium Curve Shift:
	anv = 53;
	bnv = 63;
else
	% Original:
	anv = 55;
	bnv = 65;	
end

% Pre-allocate:
if ComputationalOptions.shouldComputePeaks
	spikeFrequency_last = zeros(num_iterations, 1);
	spikeFrequency_mean = zeros(num_iterations, 1);
	spikeCounts = zeros(num_iterations, 1);
end

% For-loop form:
for i = 1:num_iterations

	% set applied current pulse
	basei = 0.0;   % baseline applied current
	pulsei = pulsei_var(i);  % magnitude of current pulse above baseline (Step_up) % 6.5 was firing start point
	t_on = 100;    % time that current pulse turns on
	t_off = 1000;  % time that current pulse turns off
	
	% set up for simulation
	options = odeset('MaxStep',1);
	%HHEquations_ftn = @(t,vars)HodHuxEquations(t, vars, basei, pulsei, t_on, t_off);
	HHEquations_ftn = @(t,vars)MultiHodHuxEquations(t, vars, basei, pulsei, t_on, t_off, anv, bnv, HH_Parameters);

	% Call the ODE solver ode15s to numerically simulate
	[t,vars] = ode15s(HHEquations_ftn, tspan, ICs, options);

	% output: t = time stamp vector of size [j by 1]. vars = [j by 4] matrix with
	% column 1 = voltage Vm, column 2 = m, column 3 = h, column 4 = n
	Vm = vars(:,1);
	m = vars(:,2);
	h = vars(:,3);
	n = vars(:,4);

	if ComputationalOptions.shouldComputePeaks
		% determine spike times and interspike intervals
		[peaks, indxs]=findpeaks(Vm,'MINPEAKHEIGHT',-10);
		spiketimes=t(indxs);
		spikeintervals=diff(spiketimes);
		spikeCounts(i) = length(spiketimes);

		% Compute Frequency:
		% IPI: Inter-peak interval: the duration (in [ms]) between the peaks times.
		% ISI: Inter-spike interval: the duration (in [ms]) between the spike times.
		if ~isempty(spikeintervals)
			% Using Two Computation Styles:
			% 1) Dr. Booth uses the last IPI to define the frequency, so I've added this as an alternative frequency metric.
			last_IPI_seconds = spikeintervals(end) / 1000; % Divide by 1000 to convert from [ms] to [sec]
			spikeFrequency_last(i) = 1 ./ last_IPI_seconds;

			% 2) This is my default definintion for frequency: the average interval between all peaks/spikes.
			mean_ISI_seconds = mean(spikeintervals) / 1000; % Divide by 1000 to convert from [ms] to [sec]
			spikeFrequency_mean(i) = 1 ./ mean_ISI_seconds;
		else
			spikeFrequency_last(i) = NaN;
			spikeFrequency_mean(i) = NaN;
		end
	end
	
	% Save the time and voltage curves for each stimulation current in case we want to plot them.
% 	time_t{i} = t;
% 	voltageTraces{i} = Vm;
	
% Compute Currents:
% 	I_Na = -HH_Parameters.gNa*(m^3)*h*(Vm-HH_Parameters.eNa);
% 	I_K = -HH_Parameters.gK*(n^4)*(Vm-HH_Parameters.eK);
% 	I_Leak = -HH_Parameters.gl*(Vm-HH_Parameters.el);
 		
	outputTraces{i}.time_t = t;
	outputTraces{i}.Vm = Vm;
	outputTraces{i}.m = m;
	outputTraces{i}.h = h;
	outputTraces{i}.n = n;
	% Compute Currents:
	outputTraces{i}.I_Na = -HH_Parameters.gNa * (m.^3) .* h .* (Vm-HH_Parameters.eNa);
	outputTraces{i}.I_K = -HH_Parameters.gK * (n.^4) .* (Vm-HH_Parameters.eK);
	outputTraces{i}.I_Leak = -HH_Parameters.gl .* (Vm-HH_Parameters.el);
	
end

% Build a table from the results for easy browsing of the different variables as a function of iteration and current.
if ComputationalOptions.shouldComputePeaks	
	resultsTable = table(pulsei_var', spikeCounts, spikeFrequency_last, spikeFrequency_mean,'VariableNames',{'AppliedCurrent','spikeCounts','spikeFrequency_last','spikeFrequency_mean'});
else
	resultsTable = table(pulsei_var','VariableNames',{'AppliedCurrent'});
end
