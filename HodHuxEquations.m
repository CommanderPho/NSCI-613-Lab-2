%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define a function containing the HH equations. Input parameter 
    % pulsei that is amplitude of applied current pulse 
    function dvdt = HHEquations(t,vars,basei,pulsei,t_on,t_off)
        
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
        
        % HH equations
        dVmdt = (-gna*(m^3)*h*(Vm-ena)-gk*(n^4)*(Vm-ek)...
                -gl*(Vm-el) + basei + pulsei*heavyside(t-t_on)*heavyside(t_off-t))/c;
        dmdt = alpham(Vm)*(1-m)-betam(Vm)*m;
        dhdt = alphah(Vm)*(1-h)-betah(Vm)*h;
        dndt = alphan(Vm)*(1-n)-betan(Vm)*n;
        
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
        
        function an=alphan(x)
            anv=55;
% 			anv=53;

            an=-0.01*(x+anv)/(exp(-(x+anv)/10)-1);
        end
        
        function bn=betan(x)
            bnv=65;
%             bnv=63;
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