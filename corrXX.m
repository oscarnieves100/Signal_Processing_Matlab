function [output_ACF, tau] = corrXX(t_signal, time, mode, plots1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Autoccorelation function of a single time
% signal x(t) by assuming ergodicity and computing the time average of the
% product x(t)x(t+tau), where tau represents the lag-time. The outputs
% include the ACF itself, and a time array tau for the lag.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% t_signal = the time signal x(t) whose autocorrelation function we are
% interested in.
%
% time = either a time array for each t in the signal x(t), or the total
% time period (duration) of the signal.
%
% mode = 1 (one-sided ACF) or 2 (two-sided ACF)
%
% plots1 = 0 (no plots) or 1 (produce plots)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% output_ACF = array of values for the autocorrelation function of t_signal
%
% tau = an array of values of the lag tau
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developed by: Oscar A. Nieves
% Institution: University of Technology, Sydney.
% Last update: 07/11/2019, 3:40pm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up default values for inputs
if nargin == 3
    plots1 = 1;
    mode = 1;
elseif nargin == 2
    plots1 = 1;
    mode = 1;
elseif nargin == 1
    plots1 = 1;
    mode = 1;
    time = linspace(0,2*pi,1000);
elseif nargin == 0
    plots1 = 1;
    mode = 2;
    time = linspace(0,2*pi,1000);
    t_signal = sin(pi*time + 0.1*randn(1,length(time)));
end

Nt = length(t_signal);
if length(time) == 1
    time = linspace(0,time,1000);
end
T = max(time);
x = t_signal;
dt = abs(time(2)-time(1));

% Calculate autocorrelation
switch mode
    case 1
        output_ACF = zeros(1,Nt);
        tau = linspace(0,T,Nt);
        xx1 = x.*conj(x);
        output_ACF(1) = dt*( (xx1(1)+xx1(end))/2 + sum(xx1(2:end-1)) );
        xshift = zeros(1,Nt);
        xshift(2:end) = x(1:end-1);
        for ii = 2:Nt
            xshift(1) = 0;
            xshift(2:end) = xshift(1:end-1);
            xx = x.*conj(xshift);
            output_ACF(ii) = dt*( (xx(ii)+xx(end))/2 + sum(xx(ii:end-1)) );
        end
    case 2
        output_ACF = zeros(1,2*Nt+1);
        tau = linspace(-T,T,2*Nt+1);
        xx1 = x.*conj(x);
        md = median(linspace(1,2*Nt+1,2*Nt+1));
        output_ACF(md) = dt*( (xx1(1)+xx1(end))/2 + sum(xx1(2:end-1)) );
        
        % Right-shift
        xshift = zeros(1,Nt);
        xshift(2:end) = x(1:end-1);
        for ii = 2:Nt
            xshift(1) = 0;
            xshift(2:end) = xshift(1:end-1);
            xx = x.*conj(xshift);
            output_ACF(md+ii-1) = dt*( (xx(ii)+xx(end))/2 + sum(xx(ii:end-1)) );
        end
        
        % Left-shift
        xshift = zeros(1,Nt);
        xshift(1:end-1) = x(2:end);
        for ii = 2:Nt
            xshift(end) = 0;
            xshift(1:end-1) = xshift(2:end);
            xx = x.*conj(xshift);
            output_ACF(md-ii+1) = dt*( (xx(1)+xx(end-ii+1))/2 + sum(xx(2:end-ii+1)) );
        end
end

switch plots1
    case 0
        disp('no plots');
    case 1
        figure(1);
        set(gcf,'color','w');
        subplot(211);
        plot(time,x,'k','LineWidth',3);
        title('Original Signal x(t)');
        xlabel('t (s)');
        ylabel('x(t)');
        xlim([min(time) max(time)]);
        ylim([min(x) max(x)]);
        set(gca,'FontSize',18);
        
        subplot(212);
        plot(tau,output_ACF,'b','LineWidth',3);
        title('Autocorrelation function R_{xx}(\tau)');
        xlabel('\tau (s)');
        ylabel('R_{xx}(\tau)');
        xlim([min(tau) max(tau)]);
        ylim([min(output_ACF) max(output_ACF)]);
        set(gca,'FontSize',18);
end
end