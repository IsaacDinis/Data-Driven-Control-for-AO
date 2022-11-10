function Kdd =  ao_dd_controller(fs,w,order,W1,W3,solver)
    %% SYS model
    Ts = 1/fs;
    
    G = tf([1],[1,0,0],Ts); % 2 samples delay
    
    %% Initial controller
    g = 0.1;
    K0 = tf([g,0],[1,-1],Ts);
    S0 = feedback(1,G*K0);
    disp(['Eigenvalues closed-loop using initial controller: ', ...
        num2str(max(abs(eig(S0)))), ' (stable CL)']) 
                                                  % < 1 --> Stable Closed-loop
    
    %% Controller order
    Fy = [1]; % fixed parts in denominator (integrator).
    Fx = 1;
    [num, den] = tfdata(K0, 'v');
    den(order + 1) = 0; % zero padding
    num(order + 1) = 0; % zero padding
    num_new = deconv(num, Fx);
    den_new = deconv(den, Fy);
    
    %% SET-UP system info
    
    [SYS, OBJ, CON, PAR] = utils.emptyStruct(); % load empty structure
    
    ctrl = struct('num',num_new,'den',den_new,'Ts',Ts,'Fx',Fx,'Fy',Fy); % assemble controller
    
    SYS.controller = ctrl;
    SYS.model = G; % Specify model(s)
    
    w_logspace = utils.logspace2(w(2),w(end),500);
    SYS.W = w_logspace;  
    
    %% Objectives

    OBJ.o2.W1 = W1;
    OBJ.o2.W3 = W3;
    
    %% Solve problem
    PAR.tol = 1e-4;
    PAR.maxIter = inf;
    % PAR.radius = 0.99;
%     tic
    [controller,obj] = datadriven(SYS,OBJ,CON,PAR,solver);
%     toc
    
    
    %% Analysis using optimal controller
    
    Kdd = utils.toTF(controller); % 
    disp(['Eigenvalues close-loop using initial controller: ', ...
        num2str(max(abs(eig(feedback(1,G*Kdd))))), ' (stable CL)']) 
                                                  % < 1 --> Stable Closed-loop
    
    %% Datadriven
    U_dd = feedback(Kdd,G);
    S_dd = feedback(1,G*Kdd);
    S_int = feedback(1,G*K0);
    opt = stepDataOptions('StepAmplitude',1);
    
    figure()
    
    subplot(2,2,1)
    step(S_dd,opt)
    title('Step Response');
    
    subplot(2,2,2)
    step(U_dd,opt)
    title('Control signal');
    
    subplot(2,2,3)
    bodemag(U_dd);
    title('Sensitivity function U');
    
    subplot(2,2,4)
    bodemag(OBJ.o2.W1^-1,S_dd,S_int,{w(2),w(end)})
    legend('W1^{-1}','Datadriven', 'Integrator');
    title('Sensitivity function S');
    
    sgtitle('Datadriven controller') 

end