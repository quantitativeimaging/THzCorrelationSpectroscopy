% Time-of-flight analysis script.

flagViscometerPlots = 1;

numberLags    = tof2-tof1 + 1;
listTaus      = zeros(numberLags, 1);
listStrength  = zeros(numberLags, 1);

for lp = 1: (numberLags);
   
    I = sum(myPow(:,(tof1-1+lp):(tof1-1+lp+smoothZlength)),2);
    
    dI = I - mean(I);
    R = xcorr(dI,'coeff'); % Note 'coeff' affects signal strength analysis
    %R = xcorr(dI);
    lags = [-length(dI)+1:length(dI)-1];
    
    myRange = fitAutoCorrRange;
    myRange = 30;

    xx      = [-myRange:myRange]';
    yy      = R( (length(dI)-myRange):(length(dI)+myRange) );

        %if(lp < 64) % neglect only for data where white noise dominates *!!
    xx(myRange+1) = []; % neglect zero lag as this is likely white noise
    yy(myRange+1) = [];
    %end
    
    % CURVE FITTING
    [xData, yData] = prepareCurveData( xx, yy );
    % Set up fittype and options.
    ft = fittype( 'gauss1' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf 0];
    opts.StartPoint = [3.34758573046279 0 1.80656622497528];
    opts.Upper = [Inf Inf Inf];
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    xxx = -50:50;
    a1 = fitresult.a1;
     if(a1 > 5 || a1 < -5)
         a1 = 0;
     end
    c1 = fitresult.c1;
    
    figure(11)
    plot(lags,R,'k','lineWidth',1)
    ylim([-0.02 1.2]);
    xlim([-100 100])
    xlabel('Lag, pulses (30 Hz)', 'fontSize', 14)
    ylabel('Autocorrelation', 'fontSize', 14)
    title('Autocorrelation of signal (XXX micron)', 'fontSize', 14)
    % title('Autocorrelation of signal (425 micron)', 'fontSize', 14)
    set(gca,'fontSize',14)

    yyy = a1*(exp(-(xxx./c1).^2));
      hold on
      plot(xxx, yyy, 'r--', 'lineWidth',2)
      %legend('Raw', 'Smoothed', 'Fitted')
      hold off
    drawnow;
    
    listTaus(lp) = fitresult.c1;
    listStrength(lp) = a1;
end


% OUTPUT: Plot the inferred velocity versus Newtonian Taylor-Couette flow
 nr = 1.4;
  %  r  = 7.5 - depths; % Assuming corrected
%  rA = 7.5;
%  rB = 3.75;
%  Ua = 9.1;   % e.g. 7.8 or 5.2  mm/s
%  vC = Ua.*( r./rA - rA./r )./(rB/rA - rA/rB);
 
 rA = 4;%3.75;% 3.75; % ANGULAR VELOCITY FORM % e.g. 3.75 mm
 % rB = 4 + (numberLags*mmPerPix)/nr; % e.g. 7 mm approximate...
 rB = 7.5;
 % Wa = 3.14;   % e.g. 1 rad/s - supplied by script that calls this one.
 gap = 3.5;
 spacing = gap/(tof2 - tof1);
 r  = rB:-spacing:rB-gap; % radial position

 vC = (-rA^2*Wa.*r + rB^2*rA^2*Wa./r)./(rB^2 - rA^2); % T.E. Faber p. 224

 % Plot velocity for Couette viscometer
 figure(10)
 % depths     = [1:numberLags]*mmPerPix;
  depths = (rB - r)*1.4; % Refractive index 1.4 --> apparent depth
 velocities = spotWidth./(listTaus / pulseRate) ;
 scatter(depths, velocities )
  xlabel('Time of Flight, mm', 'fontSize', 14)
  ylabel('Velocity, mm/s ', 'fontSize', 14)
  title('Flow Velocity', 'fontSize',14)
  set(gca,'fontSize',14)
  
 hold on

  plot(depths, vC, 'r')
 hold off
 legend('TCV measurement', 'Applied flow')