%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define a function containing the HH equations. Input parameter 
    % pulsei that is amplitude of applied current pulse 
    function dvdt = MultiHodHuxEquations(t,vars,basei,pulsei,t_on,t_off,c_anv,c_bnv, HH_Parameters)
        
        % variables:
        Vm = vars(1);
        m = vars(2);
        h = vars(3);
        n = vars(4);
        
        % set HH parameter values
		if ~exist('HH_Parameters','var')
			HH_Parameters.c=1;
			HH_Parameters.gNa=120;
			HH_Parameters.gK=36;
			HH_Parameters.gl=0.3;
			HH_Parameters.eNa=50;
			HH_Parameters.eK=-77;
			HH_Parameters.el=-54.4;
		end

		% Unpack the parameters to make things easier:
		c=HH_Parameters.c;
        gNa=HH_Parameters.gNa;
        gK=HH_Parameters.gK;
        gl=HH_Parameters.gl;
        eNa=HH_Parameters.eNa;
        eK=HH_Parameters.eK;
        el=HH_Parameters.el;
		
		% Values:
		if ~exist('c_anv','var')
			c_anv=55;
		end

		if ~exist('c_bnv','var')
			c_bnv=65;
		end

		% Additional Current Equations:
% 		I_Na = -gNa*(m^3)*h*(Vm-eNa);
% 		I_K = -gK*(n^4)*(Vm-eK);
% 		I_Leak = -gl*(Vm-el);
		
		
        % HH equations
        dVmdt = (-gNa*(m^3)*h*(Vm-eNa)... % Na+ Current
			-gK*(n^4)*(Vm-eK)... % K+ Current
            -gl*(Vm-el)... % Leak Current
			+ basei + pulsei*heavyside(t-t_on)*heavyside(t_off-t))/c; % Applied Current;
        dmdt = alpham(Vm)*(1-m)-betam(Vm)*m;
        dhdt = alphah(Vm)*(1-h)-betah(Vm)*h;
        dndt = alphan(Vm, c_anv)*(1-n)-betan(Vm, c_bnv)*n;
        
        % set output vector
        dvdt=[dVmdt; dmdt; dhdt; dndt];
    
        function am=alpham(x)
            am=-0.1*(x+40)/(exp(-(x+40)/10)-1);
        end
        
        function bm=betam(x)
            bm=4*exp(-(x+65)/18);
        end
        
        function ah=alphah(x)
            ah=0.07*exp(-(x+65)/20);
        end
        
        function bh=betah(x)
            bh=1.0/(exp(-(x+35)/10)+1);
        end
        
        function an=alphan(x, anv)
            % anv=55;
            an=-0.01*(x+anv)/(exp(-(x+anv)/10)-1);
        end
        
        function bn=betan(x, bnv)
            % bnv=65;
            bn=0.125*exp(-(x+bnv)/80);
        end
        
        function hside=heavyside(x)
            if x >= 0
                hside = 1;
            else
                hside = 0;
            end
        end
        
    end