function Kdd =  ao_dd_controller(fs,w,order,W1,W3,W32,solver)
    %% SYS model
    Ts = 1/fs;
    
    G = tf([1],[1,0,0],Ts); % 2 samples delay
    
    %% Initial controller
    g = 0.25;
    K0 = tf([g,0],[1,-1],Ts);
    S0 = feedback(1,G*K0);
    disp(['Eigenvalues closed-loop using initial controller: ', ...
        num2str(max(abs(eig(S0)))), ' (stable CL)']) 
                                                  % < 1 --> Stable Closed-loop
    
    %% Controller order
    Fy = [1,-1]; % fixed parts in denominator (integrator).
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
    
    % w_logspace = utils.logspace2(w(2),w(end),200);
    % w_logspace = linspace(w(5),w(end),200);
    w_logspace = linspace(w(1),w(end),200);
    SYS.W = w_logspace;  
    
    %% Objectives

    OBJ.o2.W1 = W1;
    OBJ.oinf.W2 = W3;
%     OBJ.o2.W3 = W32;
%     CON.W3 = tf(0.99,1,1/fs);
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
    U_int = feedback(K0,G);
    S_dd = feedback(1,G*Kdd);
    S_int = feedback(1,G*K0);
    opt = stepDataOptions('StepAmplitude',1);
    
    figure()
    
    subplot(2,2,1)
    step(S_dd,S_int,opt)
    title('Step Response');
    legend('Datadriven', 'Integrator');

    subplot(2,2,2)
    step(U_dd,U_int,opt)
    title('Control signal');
    legend('Datadriven', 'Integrator');

    subplot(2,2,3)
    bodemag(U_dd,U_int,{w(2),w(end)});
    title('Sensitivity function U');
    legend('Datadriven', 'Integrator');

    subplot(2,2,4)
    bodemag(S_dd,S_int,W1^-1,{w(2),w(end)})
    legend('Datadriven', 'Integrator','W1^{-1}');
    title('Sensitivity function S');
    
    sgtitle('Datadriven controller') 

    %%
    figure()
    bodemag(S_dd,S_int,W1^-1,{w(2),w(end)})
    legend('Datadriven', 'Integrator','disturbance^{-1}');
    title('Sensitivity function S');
    grid on;
    sgtitle('Datadriven controller') 


end