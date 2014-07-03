%% IEC 61000-4-2 Curent Discharge Waveform for 4kV
% Code written by Sebastian Courtney
% Last updated: July 1, 2014

%% Clear out figures/variables and close open windows
clear all;
close all;
clc;

%% Enable Plots (1 = enable, 0 = disable)
firstPlot = 1;
secondPlot = 1;

%% Set time range and time step. (Modify later for user input?)
step = 0.001; % time step
tFin = 80; % end of time range
t = [0:step:tFin]; % 0<ns> to tFin<ns> in 0.001<ns> steps

%% Set equation inputs (Modify later for user input?)
n = 1.8;
I1 = 16.6; % [Amps] at 4kV
I2 = 9.3; % [Amps] at 4kV
T1 = 1.1; % Tau 1 [ns]
T2 = 2; % Tau 2 [ns]
T3 = 12; % Tau 3 [ns]
T4 = 37; % Tau 4 [ns]

%% Calculate constants k1 and k2
k1 = exp(-1*(T1/T2)*(n*(T2/T1))^(1/n)); % k1 constant
k2 = exp(-1*(T3/T4)*(n*(T4/T3))^(1/n)); % k2 constant

%% Generate waveform utilizing equation on Page 13 of IEC 61000-4-2
I = (I1.*((t./T1).^n).*exp(-t./T2))./(k1*(1+(t./T1).^n)) + (I2.*((t./T3).^n).*exp(-t./T4))./(k2*(1+(t./T3).^n));

%% Generate waveform components
Ia = I1.*((t./T1).^n).*exp(-t./T2)./(k1*(1+(t./T1).^n));
Ib = I2.*((t./T3).^n).*exp(-t./T4)./(k2*(1+(t./T3).^n));

%% Find and define first peak and its index.
[Ipk,pkIndex] = max(I);

%% Find and define indices for the points equivalent to 10% and 90% of Ipk
for k=1:pkIndex
    if I(k)<=(0.1*Ipk) && I(k+1)>(0.1*Ipk)
        tenpcnt = k; % Used for shifting the 30ns and 60ns measurements to correct "zero-reference"
    elseif I(k)<=(0.9*Ipk) && I(k+1)>(0.9*Ipk)
        ninetypcnt = k;
    end
end

tRise = t(ninetypcnt) - t(tenpcnt); % Calculate Rise Time using the previously defined points

%% Generate 2nd time range that starts at t(tenpcnt), for accurate "t=0" for 30ns and 60ns plots
t2 = [t(tenpcnt):step:(tFin+((tenpcnt-1)*step))];

difT = length(t) - length(t2);

%% Find and define index for second peak
for j=(pkIndex+1):(length(I)-1)
        if I(j)>I(j-1) && I(j)>I(j+1)
            pkIndex2 = j;
            break;
        end
end

%% Plotting Code

if firstPlot==1
% Set up figure for full screen, and to append new plots/stems
figure(1);
%set(1,'units','normalized','outerposition',[0 0 1 1]);
hold on

plot(t,I); % Plot main waveform onto figure

% Plot and label the first peak [Red]
stem(t(pkIndex),Ipk,'fill','--or');
text(t(pkIndex+(1/step)),Ipk-0.1,strcat('Ipk = ',num2str(Ipk)));

% Plot and label the second peak [Red]
stem(t(pkIndex2),I(pkIndex2),'fill','--or');
text(t(pkIndex2+(0.5/step)),I(pkIndex2)+0.2,strcat('Ipk2 = ',num2str(I(pkIndex2))));

% Plot and label the 30ns and 60ns points. [Green]
stem(t((30/step)+tenpcnt),I((30/step)+tenpcnt),'fill','--og');
text(t((30.5/step)+tenpcnt),I(30/step)+0.2,num2str(I((30/step)+tenpcnt)));
stem(t((60/step)+tenpcnt),I((60/step)+tenpcnt),'fill','--og');
text(t((60.5/step)+tenpcnt),I(60/step)+0.2,num2str(I((60/step)+tenpcnt)));

% Print Rise Time to the figure and plot point at 10% of Ipk [Green]
text(t(50/step),(Ipk*(2/3)),strcat('Rise Time =  ',num2str(tRise),' ns'));
stem(t(tenpcnt),I(tenpcnt),'fill','--og');

% Plot and label the points at t=Tau1:Tau4 (T1:T4) [Cyan]
stem(T1,I(T1/step),'fill',':om');
text(T1-0.85,I(T1/step)+0.2,'T1');
stem(T2,I(T2/step),'fill',':om');
text(T2+1,I(T2/step),'T2');
stem(T3,I(T3/step),'fill',':om');
text(T3-1.5,I(T3/step),'T3');
stem(T4,I(T4/step),'fill',':om');
text(T4,I(T4/step)+0.25,'T4');

% Set more figure properties
title('IEC61000-4-2 Contact Discharge Current Waveform');
xlabel('Time (ns)');
ylabel('Discharge Current (A)');
end

%% Second Figure - Waveform Components

if secondPlot==1
figure(2)
%set(2,'units','normalized','outerposition',[0 0 1 1]);
hold on
title('IEC61000-4-2 Contact Discharge Current Waveform Components');
xlabel('Time (ns)');
ylabel('Discharge Current (A)');

plot(t,Ia,'--r');
plot(t,Ib,'--c');
plot(t,I,'--b');

stem(T1,I(T1/step),'fill',':om');
text(T1-0.85,I(T1/step)+0.2,'T1');
stem(T2,I(T2/step),'fill',':om');
text(T2+1,I(T2/step),'T2');
stem(T3,I(T3/step),'fill',':om');
text(T3-1.5,I(T3/step),'T3');
stem(T4,I(T4/step),'fill',':om');
text(T4,I(T4/step)+0.25,'T4');

stem(t((30/step)+tenpcnt),I((30/step)+tenpcnt),'fill','--og');
text(t((30.5/step)+tenpcnt),I(30/step)+0.2,num2str(I((30/step)+tenpcnt)));
stem(t((60/step)+tenpcnt),I((60/step)+tenpcnt),'fill','--og');
text(t((60.5/step)+tenpcnt),I(60/step)+0.2,num2str(I((60/step)+tenpcnt)));

stem(t(pkIndex2),I(pkIndex2),'fill','--ob');
text(t(pkIndex2+(0.5/step)),I(pkIndex2)+0.2,strcat('Ipk2 = ',num2str(I(pkIndex2))));
end

%% Create data arrays for spectre simulation

for i=1:1:2*length(t)
    if mod(i,2)~=0
        S(i) = t((i+1)/2);
        S(i+1) = I((i+1)/2);
    end
end

spectre = zeros(length(t),2);
for i=1:length(t)
    %spectre(i,1) = '+';
    spectre(i,1) = t(i);
    spectre(i,2) = I(i);
end
