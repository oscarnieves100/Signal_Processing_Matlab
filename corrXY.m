function [output_corr, tau] = corrXY(t_signal1, t_signal2, time, mode, plots1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the correlation function of a two signals x(t) and
% xy(t) by assuming ergodicity and computing the time average of the
% product x(t)y(t+tau), where tau represents the lag-time. The outputs
% include the ACF itself, and a time array tau for the lag.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% t_signal1 = the time signal x(t)
%
% t_signal2 = the time signal y(t)
%
% time = either a time array for each t in the signal x(t), or the total
% time period (duration) of the signal.
%
% mode = 1 (one-sided ACF) or 2 (two-sided ACF)
%
% plots1 = 0 (no plots) or 1 (produce plots)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% output_corr = array of values for the correlation between x(t) and y(t)
%
% tau = an array of values of the lag tau
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developed by: Oscar A. Nieves
% Institution: University of Technology, Sydney.
% Last update: 07/11/2019, 3:40pm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up default values for inputs
if nargin == 4
    plots1 = 0;
elseif nargin == 3
    plots1 = 1;
    mode = 1;
elseif nargin == 2
    plots1 = 1;
    mode = 1;
    time = linspace(0,2*pi,1000);
elseif nargin == 1
    plots1 = 1;
    mode = 1;
    time = linspace(0,2*pi,1000);
    t_signal2 = sin(pi*time + 0.1*randn(1,length(time)));
elseif nargin == 0
    plots1 = 1;
    mode = 2;
    time = linspace(0,2*pi,1000);
    t_signal2 = sin(pi*time + 0.1*randn(1,length(time)));
    t_signal1 = sin(2*pi*time + 0.3*randn(1,length(time)));
end

Nt = length(t_signal1);
if length(time) == 1
    time = linspace(0,time,1000);
end
T = max(time);
x = t_signal1;
y = t_signal2;
dt = abs(time(2)-time(1));

% Calculate correlation
switch mode
    case 1
        output_corr = zeros(1,Nt);
        tau = linspace(0,T,Nt);
        xx1 = x.*conj(y);
        output_corr(1) = dt*( (xx1(1)+xx1(end))/2 + sum(xx1(2:end-1)) );
        yshift = zeros(1,Nt);
        yshift(2:end) = y(1:end-1);
        for ii = 2:Nt
            yshift(1) = 0;
            yshift(2:end) = yshift(1:end-1);
            xx = x.*conj(yshift);
            output_corr(ii) = dt*( (xx(ii)+xx(end))/2 + sum(xx(ii:end-1)) );
        end
    case 2
        output_corr = zeros(1,2*Nt+1);
        tau = linspace(-T,T,2*Nt+1);
        xx1 = x.*conj(y);
        md = median(linspace(1,2*Nt+1,2*Nt+1));
        output_corr(md) = dt*( (xx1(1)+xx1(end))/2 + sum(xx1(2:end-1)) );
        
        % Right-shift
        yshift = zeros(1,Nt);
        yshift(2:end) = y(1:end-1);
        for ii = 2:Nt
            yshift(1) = 0;
            yshift(2:end) = yshift(1:end-1);
            xx = x.*conj(yshift);
            output_corr(md+ii-1) = dt*( (xx(ii)+xx(end))/2 + sum(xx(ii:end-1)) );
        end
        
        % Left-shift
        yshift = zeros(1,Nt);
        yshift(1:end-1) = y(2:end);
        for ii = 2:Nt
            yshift(end) = 0;
            yshift(1:end-1) = yshift(2:end);
            xx = x.*conj(yshift);
            output_corr(md-ii+1) = dt*( (xx(1)+xx(end-ii+1))/2 + sum(xx(2:end-ii+1)) );
        end
end

switch plots1
    case 0
        disp('no plots');
    case 1
        MaxV = unique( max([x y]) );
        MinV =  unique( min([x y]) );
        figure(1);
        set(gcf,'color','w');
        subplot(211);
        plot(time,x,'k',time,y,'r','LineWidth',2);
        legend('x(t)','y(t)');
        title('Original Signals');
        xlabel('t (s)');
        ylabel('Signal Amplitude');
        xlim([min(time) max(time)]);
        ylim([MinV MaxV]);
        set(gca,'FontSize',18);
        
        subplot(212);
        plot(tau,output_corr,'b','LineWidth',3);
        title('Correlation function R_{xy}(\tau)');
        xlabel('\tau (s)');
        ylabel('R_{xy}(\tau)');
        xlim([min(tau) max(tau)]);
        ylim([min(output_corr) max(output_corr)]);
        set(gca,'FontSize',18);
        clc;
end
end