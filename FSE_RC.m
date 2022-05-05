clc
clear all

a = 9; % amplitude of square signal
b = 2; % end point of first cycle square
c = 1; % sum of this and 'b' gives the period of the square signal

%THIS IS JUST A TEST OF THE FUNCTION DEFINED BELOW
%GENERATING THE SQUARE WAVE TO TO BE REPRESENTED WITH THE FOURIER SERIES EXPANSION(FSE)
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

N = 50; %summation end limit
T = b+c; %period of the signal
f = real_coef_fourier_expansion(N,t,T,s)

figure()
plot(t,s,"k")
hold on
plot(t,f)

%CREATING THE FSE WITH REAL COEFFICIENTS ALGORITHM
function f = real_coef_fourier_expansion(N,t,T,signal)
t0 = linspace(0,T,numel(t)); %integration limit
a0 = ((2/T)*trapz(t0,-signal))/2; %calculating value of a0 using the trapz function
f = a0; %assigning a0 to first element of the fourier series expansion
for k=0:N %summation with end limit = N
    fc = signal.*cos(2*pi*(1/T)*k*t); %multiply the signal with the cosine function
    fs = signal.*sin(2*pi*(1/T)*k*t); %multiply the signal with the sine function

    ak = (2/T) * trapz(t0,fc); %calculating value of ak using the trapz function
    bk = (2/T) * trapz(t0,fs); %calculating value of bk using the trapz function
    f = f+ak*cos(2*pi*(1/T)*k*t) + bk*sin(2*pi*(1/T)*k*t); %assigning the values of ak and bk to the fourier expansion
end
end
