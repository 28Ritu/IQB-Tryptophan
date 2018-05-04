% mRNA production levels as a function of intracellular tryptophan level
function [t,y] = tryptophan2()
    tspan = [0 20000];                                                % set time interval
    y1_0 = 100;                                                     % set mRNA initial concentration (in nM)
    y2_0 = 10;                                                      % set low initial tryptophan concentration = 10 nM 
    y3_0 = 200;                                                     % set trypR initial concentration (in nM)
    y4_0 = 0;                                                       % set tryp-trypR initial concentration (in nM)
    [t,y] = ode15s( @tryptophan2 ,tspan ,[y1_0, y2_0, y3_0, y4_0]);  %tryptophan evaluates r.h.s. of the ode
    
    subplot(2, 2, 1);
    plot(t, y(:,1), '-o');
    xlabel('Time');
    ylabel('mRNA Production Levels');
    title('Low Tryptophan Concentration');
    subplot(2, 2, 3);
    plot(y(:,2), y(:,1), '-o');
    xlabel('Tryptophan Production levels');
    ylabel('mRNA Production Levels');
    title('Low Tryptophan Concentration');
    
    S = 'mRNA Production Levels When Tryptophan Concentration Is Low (=10 nM)';
    disp(S);
    S1 = '      Time  mRNA Concentration (in nM)';
    disp(S1);
    disp([t, y(:,1)]);                                              % displays t and y(t)
    
    y2_0 = 500;                                                     % set low initial tryptophan concentration = 10 nM 
    [t,y] = ode15s( @tryptophan2 ,tspan ,[y1_0, y2_0, y3_0, y4_0]);  %tryptophan evaluates r.h.s. of the ode
    
    subplot(2, 2, 2);
    plot(t, y(:,1), '-o');
    xlabel('Time');
    ylabel('mRNA Production Levels');
    title('High Tryptophan Concentration');
    subplot(2, 2, 4);
    plot(y(:,2), y(:,1), '-o');
    xlabel('Tryptophan Production Levels');
    ylabel('mRNA Production Levels');
    title('High Tryptophan Concentration');
    
    S = 'mRNA Production Levels When Tryptophan Concentration Is High (=500 nM)';
    disp(S);
    S1 = '      Time  mRNA Concentration (in nM)';
    disp(S1);
    disp([t, y(:,1)]);
        function dydt = tryptophan2(t,y)
            dydt = zeros(4,1);                                                  % a column vector
            Ks = 1e-2;                                                         % gene transcription rate
            Kd = 2e-4;                                                         % mRNA degradation rate
            alpha = 1e-3;                                                       % tryp production effective kinetic constant
            Kloss = 5e-2;                                                      % tryptophan loss (protein production / degradation ....)
            Kon = 1e-4;                                                        % tryp-trypR binding
            Koff = 1e-2;                                                       % tryp-trypR dissociation
            Ng = 2;
            
            dydt(1) = Ks*Ng*(1-y(4)/200) - Kd*y(1);                           % y(1) --> mRNA
            dydt(2) = alpha*y(1) - Kloss*y(2) - Kon*y(2)*y(3) + Koff*y(4);   % y(2) --> tryptophan (tryp)
            dydt(3) = -Kon*y(2)*y(3) + Koff*y(4);                             % y(3) --> trypR (repressor)
            dydt(4) = Kon*y(2)*y(3) - Koff*y(4);                              % y(4) --> tryp-trypR complex
        end
end
