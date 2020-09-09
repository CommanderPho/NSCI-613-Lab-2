%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define a function containing the HH equations. Input parameter 
    % pulsei that is amplitude of applied current pulse 
    function dvdt = MultiHodHuxEquations(t,vars,basei,pulsei,t_on,t_off,c_anv,c_bnv)
        
        % variables:
        Vm = vars(1);
        m = vars(2);
        h = vars(3);
        n = vars(4);
        
        % set HH parameter values
        c=1;
        gna=120;
        gk=36;
        gl=0.3;
        ena=50;
        ek=-77;
        el=-54.4;

		% Values:
		if ~exist('c_anv','var')
			c_anv=55;
		end
		% 			c_anv=53;

		if ~exist('c_bnv','var')
			c_bnv=65;
        %             c_bnv=63;
		end

        % HH equations
        dVmdt = (-gna*(m^3)*h*(Vm-ena)-gk*(n^4)*(Vm-ek)...
                -gl*(Vm-el) + basei + pulsei*heavyside(t-t_on)*heavyside(t_off-t))/c;
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