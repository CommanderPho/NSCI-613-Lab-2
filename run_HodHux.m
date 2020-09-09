% script to numerically simulate the HH model
% calls function HHEquations.m which must be in same directory as this file

% set simulation timespan in ms
tspan = [0 1100];

% set initial conditions [V(t=0); m(t=0); h(t=0); n(t=0)]
v0 = -64.9997;
m0 = 0.0529;
h0 = 0.5961;
n0 = 0.3176;
ICs = [v0; m0; h0; n0];

% For-loop form:

%pulsei_var = 2:0.1:7;  % for 2(a-b): magnitude of current pulse above baseline (Step_up)
% pulsei_var = 8:1:140;  % for 3: magnitude of current pulse above baseline (Step_up)

% pulsei_var = 2:0.1:101;  % for 3: magnitude of current pulse above baseline (Step_up)

%pulsei_var = 1:1:240;  % for 4: magnitude of current pulse above baseline (Step_up)

pulsei_var = 60:1:180;  % REDO: All

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
% spikeFrequency = zeros(num_iterations, 1);
spikeFrequency_last = zeros(num_iterations, 1);
spikeFrequency_mean = zeros(num_iterations, 1);
spikeCounts = zeros(num_iterations, 1);

for i = 1:num_iterations

	% set applied current pulse
	basei = 0.0;   % baseline applied current
	pulsei = pulsei_var(i);  % magnitude of current pulse above baseline (Step_up) % 6.5 was firing start point
	t_on = 100;    % time that current pulse turns on
	t_off = 1000;  % time that current pulse turns off
	
	% set up for simulation
	options = odeset('MaxStep',1);
	%HHEquations_ftn = @(t,vars)HodHuxEquations(t, vars, basei, pulsei, t_on, t_off);
	HHEquations_ftn = @(t,vars)MultiHodHuxEquations(t, vars, basei, pulsei, t_on, t_off, anv, bnv);

	% Call the ODE solver ode15s to numerically simulate
	[t,vars] = ode15s(HHEquations_ftn, tspan, ICs, options);

	% output: t = time stamp vector of size [j by 1]. vars = [j by 4] matrix with
	% column 1 = voltage Vm, column 2 = m, column 3 = h, column 4 = n
	Vm = vars(:,1);
	m = vars(:,2);
	h = vars(:,3);
	n = vars(:,4);

	% determine spike times and interspike intervals
	[peaks, indxs]=findpeaks(Vm,'MINPEAKHEIGHT',-10);
	spiketimes=t(indxs);
	spikeintervals=diff(spiketimes);
	spikeCounts(i) = length(spiketimes);
	
	% Compute Frequency:

	% Computation Style:
	% Mean:

	% Last:
	if ~isempty(spikeintervals)
		last_IPI_seconds = spikeintervals(end) / 1000;
		spikeFrequency_last(i) = 1 ./ last_IPI_seconds;
		
		mean_ISI_seconds = mean(spikeintervals) / 1000;
		spikeFrequency_mean(i) = 1 ./ mean_ISI_seconds;
	else
		spikeFrequency_last(i) = NaN;
		spikeFrequency_mean(i) = NaN;
	end
	

	time_t{i} = t;
	voltageTraces{i} = Vm;
	
% 	fprintf('Size of Vm: %i\n', length(Vm));
	% plot Voltage vs time
% 	figure(1)
% 	plot(t,Vm,spiketimes,peaks,'*')

end

resultsTable = table(pulsei_var', spikeCounts, spikeFrequency_last, spikeFrequency_mean,'VariableNames',{'AppliedCurrent','spikeCounts','spikeFrequency_last','spikeFrequency_mean'});

