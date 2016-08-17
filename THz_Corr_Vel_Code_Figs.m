% Terahertz Correlation Spectroscopy data analysis
% Matlab scripts for analysis of time-domain terahertz sensor data to 
% produce graphs in attached paper from sample data
%
% Reference: "Terahertz correlation spectroscopy infers particle velocity 
% and rheological properties" By Eric Rees, Ke Su, Axel Zeitler
% 
% Author: Eric Rees 2016
%
% Instructions. Run this code in Matlab (tested in 2013b). It should 
% process sample data in the same folder to produce figures. 

% Figure 2 A inset 
% Time-resolved echo intensity from focal spot, for 2 bead sizes
load('ballotini_6June_004_212_30Hz'); % Load 212 micron dia. bead data
mydat = double(sam_data) -ones(size(sam_data,1),1)*mean(double(sam_data));
myPow = mydat.^2; % Use power related to E-field squared
myPow = myPow ./ (ones(size(myPow,1),1)*std(myPow)); % Convert to S/N

load('ballotini_6June_002_425_30Hz'); % Load 425 micron dia. bead data
mydat = double(sam_data) -ones(size(sam_data,1),1)*mean(double(sam_data));
myPow2 = mydat.^2; % Use power related to E-field squared
myPow2 = myPow2 ./ (ones(size(myPow2,1),1)*std(myPow2)); % Convert to S/N

myTof = 760; % This column (time of flight bin) is the in-focus data
             % each row contains a reflected signal.

figure(1)
plot( ((1:size(myPow,1)))/30 -270 ,myPow(:,myTof) +35,'r', 'lineWidth', 2)
hold on
 plot( ((1:size(myPow2,1))-2950)/ 30,myPow2(:,myTof), 'lineWidth', 1)
hold off
xlim([0 10])
   ylim([-0.1 60.5])
   xlabel('Time / seconds', 'fontSize',12)
   ylabel('THz reflection intensity', 'fontSize',12)
   % title([], 'fontSize',14')
   legend('212 \mu{}m','425 \mu{}m ', 'location', 'N')
   set(gca,'fontSize',10)
   set(gca,'yTick',[])
   set(gca,'xTick',0:2:10)
   set(gcf,'color','w')
set(gcf,'Position',[100,100,240,200]); % 540 px wide, 450 high (720/600)


% Figure 2 A main
% Two sizes of ballotini: autocorrelation of reflection, data VS fit
I212  = sum(myPow(6000:15000,(776:786)),2);; % 212 micron data. Bin to get more data. 
dI = I212 - mean(I212);
R  = xcorr(dI);

lags    = [-length(dI)+1:length(dI)-1];
myRange = 30;
xx      = [-myRange:myRange]';
yy      = R( (length(dI)-myRange):(length(dI)+myRange) );

% CURVE FITTING
% From cftool
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
c1 = fitresult.c1;
yyy = a1*(exp(-(xxx./c1).^2));

R212 = R;
R212((end+1) / 2) = [];
mylags212 = lags;
mylags212((end+1)/2) = [];
mylags212 = mylags212((end+1)/2 : (end+1)/2 + 150);
R212 = R212((end+1)/2 : (end+1)/2 + 150);
xxx212 = xxx;
yyy212 = yyy;

% 425 nm bead data - process the same way...
I425 = sum(myPow2(:,(787:796)),2); % Bin to capture more reflection data
dI = I425 - mean(I425);
R  = xcorr(dI);

lags    = [-length(dI)+1:length(dI)-1];
myRange = 20;
xx      = [-myRange:myRange]';
yy      = R( (length(dI)-myRange):(length(dI)+myRange) );

% CURVE FITTING
% From cftool
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
c1 = fitresult.c1;
yyy = a1*(exp(-(xxx./c1).^2));

R425 = R;
R425((end+1) / 2) = [];
mylags425 = lags;
mylags425((end+1)/2) = [];
mylags425 = mylags425((end+1)/2 : (end+1)/2 + 150);
R425 = R425((end+1)/2 : (end+1)/2 + 150);
xxx425 = xxx;
yyy425 = yyy;

% Plot figure of autocorrelations:
figure(2)
  scatter(mylags212/30,R212, 50,'ko','lineWidth',2)
  % ylim([-0.02 0.6]);
  xlim([0 2])
  xlabel('Time lag / seconds', 'fontSize', 18)
  ylabel('Autocorrelation', 'fontSize', 18)
  title([], 'fontSize', 14)
  % title('Autocorrelation of signal (425 micron)', 'fontSize', 14)
hold on
  plot(xxx212/30, yyy212, 'k--', 'lineWidth',2)
  legend('Data', 'Fit', 'location', 'NW')
  
  scatter(mylags212/30,R212, 50,'ro','lineWidth',2) % Over-plot in red
  plot(xxx212/30, yyy212, 'r--', 'lineWidth',2)
  scatter(mylags425/30,R425, 50,'bo','lineWidth',2)
  plot(xxx425/30, yyy425, 'b--', 'lineWidth',1)
    
hold off
  set(gca,'fontSize',16)
  set(gca,'yTick',[])
  set(gca,'xTick',0:0.5:3)
  set(gcf,'color','w')
set(gcf,'Position',[100,100,420,350]); % 540 px wide, 450 high (720/600)
mylim = ylim;
mylim(1) = -10;
ylim(mylim)
  
  
% Figure 2 B
% Ballotini velocity figure
visc = 0.100; % 101 mPa s at 25 Celsius measured. Lab at ~22 Celsius
dRho = 1760;
g = 9.8;

diameters = [212, 425, 600, 850]; % microns
radii = diameters./2;

observedVel = [ 2.16, 3.74, 8.68]; % drop-time velocity, in mm/s
observedStd = [0.2, 0.3, 1.1];     % Std dev of five drop-time velocities
tcsVel      = [ 0.45, 2.1, 3.3, 7.7 ]; % based on value at highest S:N
tcsVelStd  = [0.1, 0.4, 0.45, 1.6]; % Std of vels near highest S:N

stokesVel = radii(end).^2 *1E-9 * (2/9) * dRho * g / visc;

figure(3)
  errorbar(radii(2:4).^2 /1E6 , observedVel, observedStd, 'd', 'MarkerSize', 6, 'MarkerFaceColor', 'b', 'Linewidth', 2)
hold on
 scatter(radii.^2/1E6, tcsVel, 120, 'ok', 'lineWidth', 2);
  plot([0, radii(end).^2]/1E6, [0 stokesVel], 'lineWidth', 2)
hold off
set(gca,'fontSize',18)
%set(gca, 'FontWeight', 'bold')
xlabel('Radius^2 / mm^2 ', 'fontSize', 18)
ylabel('Velocity / mm s^{-1}', 'fontSize', 18)
legend('Drop-time','THz Correlation','Stokes', 'location', 'NW')
%set(gca,'yTick',[])
set(gca,'xTick',0:0.05:2)
set(gcf,'color','w')    
set(gcf,'Position',[100,100,420,350]); % 540 px wide, 450 high (720/600)


% Figure 3 B
% Sonogram showing time-and-depth resolved THz echo from wide-gap
% Taylor-Couette viscometer

% Load 3 Sep data, file 5 (slow-ish rotation)
clear 
load('viscometer_3Sep_5'); 
mydat = double(sam_data) -ones(size(sam_data,1),1)*mean(double(sam_data));
myPow = mydat.^2; % Use power related to E-field squared
myPow = myPow ./ (ones(size(myPow,1),1)*std(myPow)); % Convert to S/N

ti = 2200;
tf = 2500;

lag0 = 350;
lag1 = 785;

rInner = 4;
rB = 7.5;
interval = (rB - rInner)/(lag1-lag0);

myIm = myPow(ti:tf, 350:785);
myIm = fliplr(myIm);

xData = rInner:interval:rB;
yData = 0:(1/30):(tf-ti)/30;

figure(4)
surf(xData,yData,myIm)
view([0 90])
shading interp
xlim([4 7.5])
caxis([0 20])
colorbar
xlabel('Radial reflector position / mm', 'fontSize', 18)
ylabel('Time / seconds ', 'fontSize', 18)
ylabel(colorbar, 'Reflection Intensity', 'fontSize', 18);
% title('Flow Velocity', 'fontSize',14)
set(gca,'fontSize',18)
set(gcf,'color','w')
set(gcf,'Position',[100,100,720,300]); % 720 px wide, 600 high



% Figure 3 C
% Plot Graph of particle speed inferred by THz correlation spectroscopy
% against the Newtonian Taylor-Couette fluid speed at the applied rotation
% Store T-C equation-predicted velocity profiles:

% Dataset  ; viscometer rotation rate.
% 12Aug_7  ; 1.43 radian/s
% 12Aug_8  ; 2.01 radian/s
% 12Aug_9  ; 2.51 radian/s
% 12Aug_10 ; 1.03 radian/s
% 12Aug_11 ; 1.57 radian/s
% Use time-of-flight bins 460 (+/- 10) to 825 (plus/minus 20 width)
% For tests 7-11 of 12_Aug

clear
load('viscometer_12Aug_7'); % Load a dataset
Wa = 1.43; % Applied rotation rate for this dataset.

mydat = double(sam_data) -ones(size(sam_data,1),1)*mean(double(sam_data));
myPow = mydat.^2; % Assume power related to E-field squared
myPow = myPow ./ (ones(size(myPow,1),1)*std(myPow)); % Convert to S/N

tof1 = 460;            % Quickest time-of-flight to analyse (520 / 590
tof2 = 825;            % ... up to here
smoothZlength = 4;     % Boxcar smoothing width
fitAutoCorrRange = 14; % For fitting Gaussian to autocorrelation

spotWidth = 0.250*1.4; % mm  ("w" in paper: equal to 350 microns)
pulseRate = 30;        % Hz  (30 1-dimensional line scans per second)
mmPerPix  = 10/1024;   % mm

TOF;  % RUN TIME-OF-FLIGHT RESOLVED CORRELATION SPECTROSCOPY SCRIPT
eqn_t7 = vC;         % 1.43 radian/s CAPTURE theoretical result
tcv_t7 = velocities; % CAPTURE fitted result
radii  = r;

% Dataset 2
load('viscometer_12Aug_8'); % Load a dataset
Wa = 2.01; % Applied rotation rate for this dataset.
mydat = double(sam_data) -ones(size(sam_data,1),1)*mean(double(sam_data));
myPow = mydat.^2; % Assume power related to E-field squared
myPow = myPow ./ (ones(size(myPow,1),1)*std(myPow)); % Convert to S/N
TOF;  % RUN TIME-OF-FLIGHT RESOLVED CORRELATION SPECTROSCOPY SCRIPT
eqn_t8 = vC;         % 1.43 radian/s CAPTURE theoretical result
tcv_t8 = velocities; % CAPTURE fitted result

% Dataset 3
load('viscometer_12Aug_9'); % Load a dataset
Wa = 2.51; % Applied rotation rate for this dataset.
mydat = double(sam_data) -ones(size(sam_data,1),1)*mean(double(sam_data));
myPow = mydat.^2; % Assume power related to E-field squared
myPow = myPow ./ (ones(size(myPow,1),1)*std(myPow)); % Convert to S/N
TOF;  % RUN TIME-OF-FLIGHT RESOLVED CORRELATION SPECTROSCOPY SCRIPT
eqn_t9 = vC;         % 1.43 radian/s CAPTURE theoretical result
tcv_t9 = velocities; % CAPTURE fitted result

% Dataset 4
load('viscometer_12Aug_10'); % Load a dataset
Wa = 1.03; % Applied rotation rate for this dataset.
mydat = double(sam_data) -ones(size(sam_data,1),1)*mean(double(sam_data));
myPow = mydat.^2; % Assume power related to E-field squared
myPow = myPow ./ (ones(size(myPow,1),1)*std(myPow)); % Convert to S/N
TOF;  % RUN TIME-OF-FLIGHT RESOLVED CORRELATION SPECTROSCOPY SCRIPT
eqn_t10 = vC;         % 1.43 radian/s CAPTURE theoretical result
tcv_t10 = velocities; % CAPTURE fitted result

% Dataset 5
load('viscometer_12Aug_11'); % Load a dataset
Wa = 1.57; % Applied rotation rate for this dataset.
mydat = double(sam_data) -ones(size(sam_data,1),1)*mean(double(sam_data));
myPow = mydat.^2; % Assume power related to E-field squared
myPow = myPow ./ (ones(size(myPow,1),1)*std(myPow)); % Convert to S/N
TOF;  % RUN TIME-OF-FLIGHT RESOLVED CORRELATION SPECTROSCOPY SCRIPT
eqn_t11 = vC;         % 1.43 radian/s CAPTURE theoretical result
tcv_t11 = velocities; % CAPTURE fitted result


figure(5)
% Scatterplot THz correlation - but exclude outliers near the slow wall
% which arise from the wobbling wall (t_9, t_8, mainly). Discuss.
r9 = radii(tcv_t9<(eqn_t9'*2));
v9 = tcv_t9(tcv_t9<(eqn_t9'*2));
r7 = radii(tcv_t7<(eqn_t7'*2));
r8 = radii(tcv_t8<(eqn_t8'*2));
r10 = radii(tcv_t10<(eqn_t10'*2));
r11 = radii(tcv_t11<(eqn_t11'*2));
v7 = tcv_t7(tcv_t7<(eqn_t7'*2));
v8 = tcv_t8(tcv_t8<(eqn_t8'*2));
v10 = tcv_t10(tcv_t10<(eqn_t10'*2));
v11 = tcv_t11(tcv_t11<(eqn_t11'*2));

 scatter(r9(1:2:end), v9(1:2:end) , 'ko', 'LineWidth', 2 )
 hold on
   plot(radii, eqn_t9, 'k', 'lineWidth',2)
   scatter(r7(1:2:end), v7(1:2:end) , 'ro', 'LineWidth', 2 )
   scatter(r8(1:2:end), v8(1:2:end) , 'go', 'LineWidth', 2 )
   scatter(r10(1:2:end), v10(1:2:end) , 'bo', 'LineWidth', 2 )
   scatter(r11(1:2:end), v11(1:2:end) , 'mo', 'LineWidth', 2 )
   plot(radii, eqn_t7, 'r', 'lineWidth',2)
   plot(radii, eqn_t8, 'g', 'lineWidth',2)
   plot(radii, eqn_t10, 'b', 'lineWidth',2)
   plot(radii, eqn_t11, 'm', 'lineWidth',2)
   
 hold off
  xlabel('Radial position / mm', 'fontSize', 16)
  ylabel('Velocity / mm s^{-1} ', 'fontSize', 16)
  % title('Flow Velocity', 'fontSize',14)
  set(gca,'fontSize',14)
  legend('THz Correlation Spectroscopy','Newtonian Taylor-Couette Flow')
  
  set(gcf,'color','w')
  set(gcf,'Position',[100,100,540,450]); % 720 px wide, 600 high
  
  ylim([0 10])
  xlim([4 7.5])
  
  

% Figure 4
% Sonogram of time-and-depth resolved THz echo from rising bubbles in oil
clear
load('bubbles_12Aug_13'); % Load 212 micron dia. bead data
mydat = double(sam_data) -ones(size(sam_data,1),1)*mean(double(sam_data));
myPow = mydat.^2; % Use power related to E-field squared
myPow = myPow ./ (ones(size(myPow,1),1)*std(myPow)); % Convert to S/N

myIm2 = myPow(5700:6100, 630:760);
  
xData2 = 0:(10/1024):1300/1024;
yData2 = 0:(1/30):(400)/30;

figure(6)
surf(xData2,yData2,myIm2)
view([0 90])
shading interp 
xlim([0 1.25])
ylim([0 12])
xlabel('Reflector position / mm', 'fontSize', 18)
ylabel('Time / seconds ', 'fontSize', 18)
set(gca,'fontSize',18)
set(gcf,'color','w')
set(gcf,'Position',[100,100,700,300]); % 720 px wide, 600 high

% % Use cubehelix colormap (D.A. Green) for greyscale compatibility.
% % http://astron-soc.in/bulletin/11June/289392011.pdf
% mycbar = colorbar;
% colormap(cubehelix)
% colorbar('YTickLabel',[]) % Try manually edit: blank YTickLabel field
caxis([0 20])
ylabel(colorbar, 'Reflection Intensity', 'fontSize', 18);


% Figure 4 B 
% Graphs of a single bubble reflection 
mybub = myIm2(210:270, 1:30);

figure(7)
plot((0:(10/1024):(30-1)/102.4),sum(mybub,1), 'lineWidth', 2);

xlim([0 0.3])
ylim([0 310])
set(gca,'yTick',[])
xlabel('Reflector position / mm', 'fontSize', 18)
ylabel('Intensity ', 'fontSize', 18)
set(gca,'fontSize',16)
set(gcf,'color','w')
set(gcf,'Position',[100,100,300,240]); % 720 px wide, 600 high

figure(8)
plot((210/30:(1/30):(270)/30),smooth(sum(mybub,2),3), 'lineWidth', 2);
% xlim([0 0.3])
ylim([0 120])
set(gca,'yTick',[])
xlabel('Time / seconds', 'fontSize', 18)
ylabel('Intensity ', 'fontSize', 18)
set(gca,'fontSize',16)
set(gcf,'color','w')
set(gcf,'Position',[100,100,300,240]); % 720 px wide, 600 high


% End of Figure 4

% Figure 1C 
% Supplementary panel showing waterfall plot of successive measurements
% 
clear
load('ballotini_6June_002_425_30Hz'); % Load 425 micron dia. bead data
mydat = double(sam_data) -ones(size(sam_data,1),1)*mean(double(sam_data));
myPow2 = mydat.^2; % Use power related to E-field squared

% figure(14)
% for lp = 6210:2:6240
%    
%     plot([1:1024]/100 -5,myPow2(lp,:),'b', 'lineWidth', 2);
%     hold on
%       plot([710:764]/100 -5,myPow2(lp,710:764),'r', 'lineWidth', 2);
%       % plot([610:664]/100 -5,sam_data(lp,610:664),'k', 'lineWidth', 2);
%     hold off
%         ylim([0 0.5])
%     
%     xlabel('~ Time of flight / "mm" ', 'fontSize', 14)
%     ylabel('Reflected THz E-field', 'fontSize', 14)
%     
%     pause(1)
%     % title(['425 micron silica, ', num2str( ((lp-1710)*1000/30),4 ),...
%     %       'ms'], 'fontSize', 14)
% 
%     % myFig = getframe(gcf); 
%     % myIm = myFig.cdata;
%     % imwrite(myIm,['C:\Documents and Settings\ejr36\My Documents\', ...
%     %               'Work_Papers\2014_THz_Corr_Vel\PPT\anim\im', ...
%     %               int2str(lp),'.png'],'png');
% end

figure(15)
for lp = 6221:1:6227
    
    hold on
      plot([500:1000]/100 -5, smooth(myPow2(lp,500:1000),7) - 0.1*lp,'k', 'lineWidth', 2);
      plot([710:764]/100 -5, smooth(myPow2(lp,710:764),7)- 0.1*lp,'r', 'lineWidth', 2);
    hold off
end

set(gca,'yTick',[])
xlim([0 4.5])
ylim([-622.8 -622])
xlabel('Time of flight / mm', 'fontSize', 16)
ylabel('Reflection intensity ', 'fontSize', 16)
set(gca,'fontSize',16)
set(gcf,'color','w')
set(gcf,'Position',[100,100,280,280]); % 720 px wide, 600 high