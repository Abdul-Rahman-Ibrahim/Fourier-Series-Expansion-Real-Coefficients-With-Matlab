clc
clear

a = 9; % fifth digit of student number
b = 2; % sixth digit of student number
c = 1; % last digit of student number
Ns = [3,5,10,50]; %values of N

%GENERATING THE SQUARE WAVE
t = 0:0.01:15; %time vector with an increment of 0.01
s = zeros(size(t)); %square wave signal
for k=0:4 %generate 5 period of square signal
    for ii=1:numel(t)
        if ((t(ii)>=k*(b+c)) && (t(ii)<=k*(b+c)+b))
            s(ii) = a;
        elseif ((t(ii)>k*(b+c)+b) && (t(ii)<(k+1)*(b+c)))
            s(ii) = 0;
        end
    end
end

%GENERATING THE TRIANGULAR WAVE
t = 0:0.01:20; %time vector with an increment of 0.01
v = zeros(size(t)); %traingular wave vector
T = 4; %period
for k=0:4 %generate 5 period of triangular signal
    for ii=1:numel(t)
        if ((t(ii)>=k*T) && (t(ii)<=(0.5*T + k*T)))
            v(ii) = (1 - 4/T*(t(ii)-k*T));
        elseif ((t(ii)>(0.5*T + k*T)) && (t(ii)<(T + k*T)))
            v(ii) = ((4/T*(t(ii)-k*T)) - 3);
        end
    end
end

%GENERATING THE FSE OF THE SQUARE WAVE
TITLE = 'FSE of Square Wave with real coefficients';
XLABEL = 't';
YLABEL = 's(t)';
AXIS = [0,16,0,50];
t = 0:0.01:15;
T = b+c;
t0 = linspace(0,T,numel(t)); %integration limit
for N=Ns
    %call the FSE function
    real_coef_fourier_expansion(N,t,t0,T,TITLE,XLABEL,YLABEL,AXIS,s)
end

%GENERATING THE FSE OF THE TRIANGULAR WAVE
TITLE = 'FSE of Triangular Signal with real coefficient';
XLABEL = 't';
YLABEL = 'v(t)';
AXIS = [0,20,0,2];
t = 0:0.01:20;
T = 4;
t0 = linspace(0,T,numel(t)); %integration limit
for N=Ns
    %call the FSE function
    real_coef_fourier_expansion(N,t,t0,T,TITLE,XLABEL,YLABEL,AXIS,v)
end

%CREATING THE FSE WITH REAL COEFFICIENTS ALGORITHM
function real_coef_fourier_expansion(N,t,t0,T,TITLE,XLABEL,YLABEL,AXIS,signal)
a0 = ((2/T)*trapz(t0,-signal))/2; %calculating value of a0 using the trapz function
f = a0; %assigning a0 to first element of the fourier series expansion
for k=0:N %summation with end limit = N
    fc = signal.*cos(2*pi*(1/T)*k*t); %multiply the signal with the cosine function
    fs = signal.*sin(2*pi*(1/T)*k*t); %multiply the signal with the sine function

    ak = (2/T) * trapz(t0,fc); %calculating value of ak using the trapz function
    bk = (2/T) * trapz(t0,fs); %calculating value of bk using the trapz function
    f = f+ak*cos(2*pi*(1/T)*k*t) + bk*sin(2*pi*(1/T)*k*t); %assigning the values of ak and bk to the fourier expansion
end
%create graph for the fourier series expansion
figure
plot(t,signal,'k',LineWidth=1) % plot the original signal
axis(AXIS) %sets axis
title(TITLE) %sets title
xlabel(XLABEL) %sets xlabel
ylabel(YLABEL) %sets ylabel
grid on %sets grid
hold on
plot(t,f,'r') %plotting the fourier series expansion of the signal
legend(YLABEL, "FSE(N="+num2str(N)+")") %sets legend
end