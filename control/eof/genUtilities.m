classdef genUtilities        
    methods (Static)
        function [IM, Cp, Cm, IM4p, IM4m] = IMcalibration(ngs, tel, wfs, dm, amp, calib_type, H, Nstep)
            %
            % IM: 'slope' interaction matrix
            % Cb: not normalized interaction matrix
            % IM4: meta intensities interaction matrix, if applicable
            %
            %
            % H = calibration poke matrix. Can be 'eye', 'hadamard', or an
            % inversible matrix, like a DM KL matrix
            %
            % Nstep = Number of step for performaing the parralelized calirabtion
            % =============================================================

            % ngs = ngs.*tel*dm*wfs;
            
            % Save the initial position of the DM 
            dm0 = dm.coefsDefault*0;
            % dm0 = dm.coefs;
            IM4 = [];


            if ~exist('Nstep', 'var')
                Nstep = 1;
            end

            nPushPull = 1;

            basisType = 'modal';
            switch calib_type
                case 'step'
                    disp('Step calibration')
                    basisType = 'zonal';
                    % Standard step calibration
                    Nmode = dm.nValidActuator;
                    H = eye(Nmode);                    
                case 'hadamard'
                    % Hadamard calibration
                    disp('Hadamard calibration')
                    basisType = 'zonal';
                    H = [];
                    Nmode = dm.nValidActuator;
                    while isempty(H) % Look for smallest Hadamard matrixs
                        try
                            H  = hadamard(Nmode);
                        catch
                            Nmode = Nmode + 1;
                        end
                    end    
                case 'kl'
                    Nkl    = size(H,2);
                    Nmode  = Nkl;
                otherwise
                    disp('User poke matrix')
                    Nmode = dm.nValidActuator;
            end
            NmodePerIter = floor(Nmode/Nstep)+1;
            nSlope       = wfs.nSlope;

            camres = wfs.camera.resolution;
            if strcmp(wfs.tag, 'SHACK-HARTMANN')
                camX = camres(1); camY = camres(2);
            else
                camX = size(wfs.slopesMap,1);
                camY = size(wfs.slopesMap,2);
            end
            IM4p = zeros(camX, camY, Nmode);
            IM4m = zeros(camX, camY, Nmode);

            sp = zeros(nSlope, NmodePerIter, nPushPull);
            sm = zeros(nSlope, NmodePerIter, nPushPull);

            % Start calibration
            msg = [];
            for ical=1:Nstep
                if ~isempty(msg) 
                     fprintf(repmat('\b', 1, strlength(msg)));  % Erase previous prompt line
                end
                msg = sprintf('### Calibration progress: %3.0f / %.0f',ical,Nstep);
                fprintf(msg)

                % PUSH-PULL
                i0 = (ical-1)*NmodePerIter+1;
                i1 = min(ical*NmodePerIter, Nmode);
                modeList = i0:i1;

                Hvec = H(1:dm.nValidActuator,modeList)*amp;
                for kf = 1:nPushPull

                    % PUSH
                    % dm.coefs = dm.coefsDefault;
                    dm.coefs = dm0 + Hvec;
                    +ngs;
                    
                    sp(:,modeList,kf)  = wfs.slopes;
                    if ~strcmp(wfs.tag, 'SHACK-HARTMANN')
                        IM4p(:,:,modeList) = reshape(wfs.slopesMap, camX,camY,[]);
                    end

                    % PULL
                    dm.coefs = dm0 - Hvec;
                    +ngs; 
                    sm(:,modeList,kf)  = wfs.slopes;
                    if ~strcmp(wfs.tag, 'SHACK-HARTMANN')
                        IM4m(:,:,modeList) = reshape(wfs.slopesMap, camX,camY,[]);
                    end
                end
                drawnow

            end % end of calibration loop
            dm.coefs = dm0; % Reset DM positions
            +ngs;

            Cp = mean(sp,3);
            Cm = mean(sm,3);

            fprintf('\n')
            disp('ok')
            Cb  = 0.5*(Cp-Cm);
            switch basisType
                case 'zonal'
                    IMb = Cb*inv(H);
                    IM  = IMb(:,1:dm.nValidActuator)/amp; % Cut away zeros from the resulting IM
                otherwise
                    IM = Cb/amp;
            end
        end
  
        
        function IM = pokeCalibration(ngs, tel, wfs, dm, calAmp, pokeMatrix, Nstep)
            % Calibrate system according to the given poke matrix. It can
            % be only a partial calibration.
            %
            % Only used in convolutional model for now.
            %
            % =============================================================   
            dm0 = dm.coefs;
            

            [n,Nmode] = size(pokeMatrix);
            if n~=dm.nValidActuator
                pokeMatrix = pokeMatrix';
            end
            [n,Nmode] = size(pokeMatrix);
            
            nPushPull = 1;
            
            nSlope = wfs.nSlope;
            Cp = zeros(nSlope,Nmode);
            Cm = zeros(nSlope,Nmode);

            sp = zeros(nSlope, nPushPull);
            sm = zeros(nSlope, nPushPull);
           
            nC = floor(Nmode/Nstep);
            u  = 0;
            IM = [];
            msg = [];
            istep = 0;
            while u(end)<Nmode
                istep = istep+1;
                u = u(end)+1:min(u(end)+nC,Nmode);

                if ~isempty(msg) 
                     fprintf(repmat('\b', 1, strlength(msg)));  % Erase previous prompt line
                end
                msg = sprintf('### Calibration progress: %4d:%4d / %4d',u(1),u(end),Nmode);
                fprintf(msg)

                for kf = 1:nPushPull
                    dm.coefs = pokeMatrix(:,u)*calAmp;
                    +ngs
                    sp(:,u,kf) = wfs.slopes;

                    dm.coefs = - pokeMatrix(:,u)*calAmp;
                    +ngs
                    sm(:,u,kf) = wfs.slopes;
                end
                drawnow
            end % end of calibration loop
            IM = 0.5*(mean(sp,3) - mean(sm,3))/calAmp;

            dm.coefs = dm0; % Reset DM positions
            fprintf('\n')
        end




        %% ================================================================
        
        function atm = paranalAtmosphere(seeing, photoAtm, meanWindSpeed, zenithAngleInRad)

            % Parameters of Paranal atmosphere 9-layer model from Nelly Cerru
            nlayer         = 9; % # of layers
            L0             = 25 ; % 25 [Nelly] % Turbulence external scale (m) % 22 = Eris sims
            altitude       = [0.042 0.140 0.281 0.562 1.125 2.25 4.5 9 18]*1000; % Altitude of each layer in the atmosphere(m)
            windSpeed      = [15 13 13 9 9 15 25 40 21]  ; % Wind speed on each layers in pure frozen flow
            fractionnalR0  = [53.28 1.45 3.5 9.57 10.83 4.37 6.58 3.71 6.71]/100 ; % Cn2 profile coefficient for each layer
            windDirection  = [38 34 54 42 57 48 -102 -83 -77]*pi/180; %2*pi*rand(1,param.nlayer) ; % Wind direction on each layer, compare to the ground


%             % Guido Agapito [ERIS]
            % nlayer        = 10; % # of layers
            % L0            = 22;
            % altitude      = [110, 238, 402, 728, 1381, 2686, 5296, 9066, 12836, 16316]; %m
            % fractionnalR0 = [0.59, 0.02, 0.04, 0.06, 0.01, 0.05, 0.09, 0.04, 0.05, 0.05];
            % windSpeed     = [6.6, 5.9, 5.1, 4.5, 5.1, 8.3, 16.3, 30.2, 34.3, 17.5];
            % % windSpeed     = [6.6, 5.9, 5.1, 4.5, 5.1, 8.3, 16.3, 30.2/6, 34.3/6, 17.5/2];
            % windDirection  = [38 34 54 42 57 48 -102 -83 -77, 34]*pi/180; %2*pi*rand(1,param.nlayer) ; % Wind direction on each layer, compare to the ground

            fractionnalR0 = fractionnalR0/sum(fractionnalR0);

            % Input parameters
            if ~exist('zenithAngleInRad', 'var')
                zenithAngleInRad = 0.;
            end
            if ~exist('meanWindSpeed', 'var')
                windSpeedRatio = 1;
            else
                meanWindSpeed0 = sum(windSpeed.^(5./3).*fractionnalR0).^(3/5); % Wind speed of the sum of layers; ~15.5m/s
                windSpeedRatio = meanWindSpeed/meanWindSpeed0;
                windSpeed      = windSpeed*windSpeedRatio;
            end
            if ~exist('photoAtm', 'var')
               photoAtm = photometry.V; 
            end            
            if ~exist('seeing', 'var')
                seeing = 0.75;
            end
            
            r0 = 550e-9/seeing/4.85e-6; % photometry.V wavelength
            r0 = r0*(photoAtm.wavelength/photometry.V.wavelength)^-1.2;          
            
            atm = atmosphere(photometry.V, r0, 'L0', L0,...
                        'altitude',      altitude,...
                        'fractionnalR0', fractionnalR0,...
                        'windSpeed',     windSpeed,...
                        'windDirection', windDirection, ...
                        'zenithAngle', zenithAngleInRad);                   
        end
        
        
        %% ===============================================================
                       
        
        function drawModes(modes, modeList, id)
            % drawModes(modes, modeList, id)
            % -------------------------------------------------------------
            % Draw the modes required by the list. The number of
            % displayed modes is limited the 25 first.
            %
            % -------------------------------------------------------------
            NmodeMax = 50;
            Nmode = min(NmodeMax, numel(modeList));
            modeList = modeList(1:Nmode);
            
            Ny = ceil(sqrt(Nmode));
            Nx = max(1, ceil(Nmode/Ny));
            fig=figure(id);
            k=0;
            
            [res1, res2, nMode] = size(modes);
            mask = sum(abs(modes),3)<1e-9;
            for i=1:Nx
                for j=1:Ny
                    k=k+1;
                    if k<=Nmode
                        mode = modes(:,:,modeList(k));
                        im = max(abs(squeeze(mode(:))),[],'all', 'omitnan');
                        mode(mask) = NaN;

                                                
                        subplot(Nx,Ny,k)                        
                        %h=subaxis(Ny,Nx,j,i,.9,.9, 'S', 0, 'P',0.);                  
                        imagesc(mode)
                        colormap('gray')
                        caxis([-im,im])
                        title(sprintf('KL #%.0f', modeList(k)))
                        axis square tight off
                    end
                end
            end
        end
       
        

        function petals = makePetals(tel, petalAngle)
            % Angles of petals in radian
            if ~exist('petalAngle', 'var')
                petalAngle = 0;
            end
            Npix = tel.resolution(1);

            petals = zeros(Npix,Npix,4);
            xx=linspace(-Npix/2,Npix/2,Npix);
            [xg,yg]=meshgrid(xx,xx);
            tg = angle((xg+1j*yg)*exp(1j*petalAngle))+pi;
            
            for i = 0:3
                m1 = tg  >     i*pi/2;
                m2 = tg <= (i+1)*pi/2;
                petals(:,:,i+1)= m1.*m2;
            end
        end


        %% ================================================================
        
        
        function [photoX, photoY] = photoCenter(psf, power)
            % [photoX, photoY] = photoCenter(psf, power)
            %--------------------------------------------------------------
            % This function computes the photocenter of a PSF. In order to
            % reduce the impact of background, the photocenter can be
            % estimated on a certain power of the initial PSF. The default 
            % power is 1. The higher the power, the closer from a 'PSF max'
            % estimator. 
            %
            % INPUT
            %  - psf: the image from which we look for the photocenter
            %  - power (optional): photocenter can be estimated on
            %  psf^power instead of psf to reduce background inflluence.
            %  But it also increases the influence of hot pixels.
            %  Default value is power=1.
            %
            % OUTPUT
            %  - [photoX, photoY]: photocenter in X and Y direction (with X = 1st dim. of psf) 
            %
            %--------------------------------------------------------------
            
            if ~exist('power', 'var')
                power = 1;
            end
            
            % Francois's version
            [Nx,Ny] = size(psf);
            sum_psf = sum(sum(psf.^power));

            photoX = (1:Nx)*sum(psf.^power,2)/sum_psf;
            photoY = sum(psf.^power,1)*((1:Ny)')/sum_psf;

        end % photoCenter
        
        %% ================================================================
        % Control
        function [H, hwfs, hsys, hol] = integratorRtf(f_loop, dit, delay, gain)            
            % Single integrator controller model
            hwfs = sinc(pi*f_loop*dit).*exp(-1j*pi*f_loop*dit)/dit;
            hsys = exp(-1j*pi*f_loop*delay)./(1j*2*pi*f_loop);
            hol  = hwfs.*hsys;
            H    = 1./(gain*hol+1);
        end


        function [amp_cor, amp_n] = AO_rtf(freq, samplingTime, gain, pureDelay, K)
            % pureDelay: delay in addition to the frame one!
                        
            % Following Gendron 1994
            h_wfs   = sinc(pi*freq*samplingTime).*exp(-1i*pi*samplingTime);
            h_delay = exp(-2i*pi*pureDelay*freq)./(2i*pi*freq)/samplingTime;
        
            h_sys   = h_delay;
            h_ol    = h_wfs.*h_sys;     
            
            %
            G_ol  = gain*h_ol;
            G_sys = gain*h_sys;
        
            % CL error transfer function
            E0_int  = 1./(1+G_ol);         % NB
            amp_cor = abs(E0_int).^2 * K;
        
            % Noise transfer function
            H0    = G_sys./(1+G_ol);       % NB                
            amp_n = abs(H0).^2 * K;
        end        


        %% ================================================================
        
        
        function [Ip, rn] = psfProfile(I, pixelScale, power)
            % Compute profile of a PSF
            %
            if ~exist('pixelScale', 'var')
                pixelScale = 1;
            end
            if ~exist('power', 'var')
                power = 1;
            end
            
            N = size(I,1);
            xg = (1:N);
            yg = (1:N);
            
            [xg, yg] = meshgrid(xg, yg);
            
            % Find photocenter of the image
            [pcx, pcy] = genUtilities.photoCenter(I,power);
            
            % Define grid for polar interpolation
            zg = xg + 1j*yg;
            zg = zg - mean(zg(:));
            zg = xg + 1j*yg - (pcy+1j*pcx);
            rg = abs(zg);
            tg = angle(zg);

            rn = linspace(0, max(rg(:)), N);
            tn = linspace(-pi,pi,33);
            [rgn, tgn] = meshgrid(rn, tn);

            Ip = griddata(rg, tg, I, rgn, tgn);            
            Ip = median(Ip, 1, 'omitnan');
            Ip(1) = max(I(:)); % Fix for missing the central value!
            
            rn = rn*pixelScale;
        end


        %% ================================================================

        
        function ee = getEnsquaredEnergy(eeWidthInDiff,psf,D,nyquistSampling)
            
            % Get the otf and the strehl
            otf    = real(tools.psf2otf(psf));
            if nyquistSampling>=1
                otf    = tools.crop(otf,round(size(otf,1)/nyquistSampling));
            end
            otf    = otf/max(otf(:));
            
            % entrapped energy
            k = 0;
            ee = [];
            for eeWidth = eeWidthInDiff
                k        = k + 1;
                a        = eeWidth/D;%(wvl/D*constants.radian2arcsec)/D;
                nOtf     = length(otf);
                u        = linspace(-1,1,nOtf).*D;
                [x,y]    = meshgrid(u);
                eeFilter = a^2*sinc(x.*a).*sinc(y.*a);
                ee(k)    = trapz(u,trapz(u,otf.*eeFilter));
            end            
        end
        
        
        
        %% ================================================================

        
        function [PSD, freq, cumulativePSD] = PSDavg(time, datat, sectionLength)
            % Function [PSD, freq, cumulativePSD] = PSDavg(time, data, sectionLength)
            %--------------------------------------------------------------
            % Compute the PSD of the data set.
            % !!! Returns the squared modulus of the FFT of the signal. !!!
            %
            % INPUT
            %  time: acquisition time of the data set 
            %    /!\ time must be uniformly sampled /!\
            %  data: data set
            %  sectionLength: length of sections used for computing PSDs,
            %  before they are averaged
            %
            % OUTPUT
            %  - Signal PSD
            %  - Frequency bins of PSD
            %  - Cumulative PSD - cPSD(end) = var(datat)
            %
            % FIX ME
            %
            % -------------------------------------------------------------
            
            Ndata     = length(datat);
            Navg      = floor(Ndata/sectionLength);  % round to next smaller integer (remove last bunch of data)
            Navg      = max(Navg, 1);  % in case the section length does not allow for Navg = 2
            if Navg < 1.5
               sectionLength = Ndata; 
            end
            
            time      = time - time(1);
            N         = sectionLength;
            dt        = mean(diff(time));
            max_ntime = N*dt;
            freq      = linspace(1./max_ntime, 1/dt, sectionLength);
            
            PSDi = zeros(Navg,N);
            for i=1:Navg
                lim_inf = (i-1)*N+1;
                lim_sup = i*N;

                ddd = datat(lim_inf:lim_sup);
                ddd = ddd - mean(ddd); % Remove mean of each section
                xdft      = fft(ddd);
                PSDi(i,:) = (1/N) * abs(xdft).^2;
            end
            
            PSD = mean(PSDi, 1);
            cumulativePSD = 2*cumsum(PSD(1,1:floor(N/2)))/N;
            
            PSD = PSD(1:floor(N/2));
            freq = freq(1:floor(N/2));
        end % PSDavg
  
       
        
        %% ================================================================
        
        function ac = autocorr(x)
            % ac = autocorr(x)
            % -------------------------------------------------------------
            % Compute the auto-correlation of a signal
            %
            % INPUT
            %   - x: the signal to autocorrelate
            %
            % OUTPUT
            %  - ac: auto-correlation of the signal
            %
            % -------------------------------------------------------------

            ac = genUtilities.crosscorr(x,x);
        end
        
        %% ================================================================ 
 
        function xc = crosscorr(x,y)
            % xc = crosscorr(x,y)
            % -------------------------------------------------------------
            % Compute the cross-correlation of two signals
            %
            % INPUT
            %   - x, y: two signals with same number of samples
            %
            % OUTPUT
            %  - xc: cross-correlation of signals
            %
            % -------------------------------------------------------------

            corrLength=length(x)+length(y)-1;
            xc = real(fftshift(ifft(fft(x/length(x),corrLength).*conj(fft(y/length(y),corrLength)))));
        end
     
        %% ================================================================

        function [xc, yc] = circle(r0)
            phi=linspace(0, 2*pi, 100);
            
            xc = r0*cos(phi);
            yc = r0*sin(phi);            
        end
        
        
        %% ================================================================
        
        
        function mask = maskGen(Nx, Ny, x0, y0, r, obscuration)
            % mask = maskGen(Ny, Nx, x0, y0, r, obscuration)
            % --- Legacy NIRPS function ---
            %--------------------------------------------------------------
            % This function generate a circular mask
            %
            % input:    Nx,Ny   - size of grid in pix
            %           x0,y0   - position of the circle center [pix]
            %           r       - radius in [pix]
            % output:   mask    - array with a circular mask. Value is 1 if 
            %                    the distance from center <= r
            %--------------------------------------------------------------
            [Xg,Yg] = meshgrid((1:Nx)-x0,(1:Ny)-y0);
            
            mask  = (abs(Xg+1j*Yg)) <= r;
            
            if exist('obscuration', 'var')
                mask2  = (abs(Xg+1j*Yg)) < obscuration*r;
                mask = mask - mask2;
            end
            
        end % maskGen
        
        function mask = cirularMask(N, D, D0, obscurationRatio)    
            %
            % D:  grid diameter
            % D0: mask diameter
            %
            x       = linspace(-D/2,D/2,N);
            [xx,yy] = meshgrid(x,x);
            rr      = abs(xx+1j*yy);
            mask    = rr<=D0/2;

            if exist('obscurationRatio', 'var')
                mask2    = rr<=obscurationRatio*D0/2;
                mask = mask - mask2;
            end
        end


        
        %% ================================================================
        % Dispersion laws
%         function R = dispersionOwen(lbd, T, P, H) 
          % SOMETHING FISHY WITH FORMULA -- MISSING LBD DEPENDENCE
%             % Owens 1967 formula
%             Pw_sat = 6.11*10.^((7.5*(T-273.15))/(T+237.7-273.15));  % pression de saturation de la vapeur d'eau
%             Pw = Pw_sat * H;    % pression d'eau d'apres l'hygrometrie
%             Ps = P - Pw;        % pression d'air
% 
%             Ds = Ps/T*(1 + Ps*(57.9e-8 - 9.325e-4/T + 0.25844/T.^2) );
%             Dw = Pw/T*(1 + Pw*(1 + 3.7e-4*Pw)*(-2.37321e-3 + 2.23366/T - 710.792/T.^2 + 7.75141e4/T.^3));
% 
%             R = (2371.34 + 683939.7/(130.-sig.^2) + 4547.3/(38.9-sig.^2))*Ds + (6487.31 + 58.058*sig.^2 - 0.7115*sig.^4 + 0.8851*sig.^6)*Dw;
%             R = R*1.e-8;
%         end

       
        function R = EDLEN(lbd, T, P)
            % # Edlen 1966 formula
            % lbd in [nm]
            R = 8342.13 + 2406030./(130. - (lbd/1000).^-2) + 15997./(38.9 - (lbd/1000).^-2); % Edlen
            R = R*0.00138823.*(P/1.33322)./(1 + 0.003671.*(T-273.15));   % to take into account the temperature and pressure
            R = R*1.e-8;
        end
        
        function dispInArcsec = dispersionAngle(T, P, zenithInDeg, lbdSci, lbdNgs)
            % T in Kelvin
            % P(ressure) in mBar
            % zenithInDeg in deg
            % lbdSci, lbdNgs in meters
            %            
            Rsci = genUtilities.EDLEN(lbdSci*1e9,T,P);
            Rngs = genUtilities.EDLEN(lbdNgs*1e9,T,P);
            
            dispInRad = tan(zenithInDeg/180.*pi).*(Rsci-Rngs);      
            dispInArcsec = dispInRad/4.85e-6;
        end
        
        function dispInArcsec = dispInArcSec_OOMAO( zenithAngleInDeg, lbdSci, lbdNgs)
            n     = @(x) (8.34213e-5 + 0.0240603./(130-x.^2)+ 0.00015997./(38.9-x.^2)); % Hardy Eq. (3.16). Index of refraction (and not refractivity) variations at standard temperature and pressure 
            dispInArcsec = (n(lbdSci*1e6) - n(lbdNgs*1e6)).*tan(zenithAngleInDeg*pi/180)/4.85e-6;
        end




        
        
        
    end % end methods

    
    % =====================================================================
    % Folder utilities
    % =====================================================================
    methods (Static)
              
        function createFolder(folder, verbose)
            % createFolder(folder)
            % -------------------------------------------------------------
            % Create a new folder. Maybe.
            %
            % INPUT
            %  - folder: the folder name (global path, or local)
            %  - verbose: If ==1 (default), it display a message if folder already exists. 
            % OUTPUT
            %  - NONE
            %
            % -------------------------------------------------------------
            if ~exist('verbose', 'var')
                verbose = 0;
            end
            
            if ~exist(folder, 'dir')
                mkdir(folder);
                disp(strcat('Folder created:    ', folder))
            else
                if verbose
                   disp(strcat(folder, ' already exists.'))
                end
            end
        end % createFolder

        %==================================================================

        function timeStamp = createTimeStamp()
            % createTimeStamp()
            % -------------------------------------------------------------
            % Create a time stamp.
            %
            % INPUT
            %  - NONE
            %
            % OUTPUT
            %  - timeStamp: string with time stamp corresponding to the time
            %    the fucntion is called
            %
            % -------------------------------------------------------------

            ts = clock();
            timeStamp = strcat(datestr(ts,'yyyy_mm_ddTHH_MM_SS'));
            
        end
         

        
    end
end
