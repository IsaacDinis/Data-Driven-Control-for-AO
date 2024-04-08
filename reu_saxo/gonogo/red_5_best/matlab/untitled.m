if sol.nIter==0
    Ts = system.controller.Ts;
    z = tf('z',Ts);

    Fx = utils.resp(tf(system.controller.Fx,1,Ts),system.W); % fixed part Y
    Fy = utils.resp(tf(system.controller.Fy,1,Ts),system.W); % fixed part Y
    
    z_ = utils.resp(z,system.W);
    Zy = z_.^((szy-1):-1:0); % [0,z,...,z^nx]
    Zx = z_.^((szx-1):-1:0); % [0,z,...,z^ny]
    ZFy = Zy.*Fy;
    ZFx = Zx.*Fx; % Frequency response K  with the fixed parts
    
    PLANT = zeros(nCon,nMod);
    for mod =  1: nMod
        PLANT(:,mod) = utils.resp(system.model(:,:,mod),system.W); % Model Frequency response
    end
    