clc
clear

% PROBLEM 2

n = -12:12; %integer vector

%signal x[n]
x = (n+4) .* ((n>=-5) & (n<=-1)) + ...
    (3) .* ((n>-1) & (n<2)) + ...
    (7-2*n) .* ((n>=2) & (n<=5)) + ...
    (0) .* (n>5);

for ii=1:4 % generate 4 graphs for each 4 different value of 'a'.
    a = get_a_value(ii); %get the current value of 'a' to plot

    %impulse response h[n]
    h = (1) .* ((a >=0) & ((n>=0) & (n<=a))) + ...
        (1) .* ((a <= 0) & ((n>=a) & (n<=0))) + ...
        (0) .* (n>0);
    
    % get convolution sum from the function defined below.
    [y,l] = linearconvolve(x,h); %get the convolution sum and the number of elements in y
    
    %generate plot for current value of 'a'
    figure()
    n = linspace(-12,12,l); %generate vector n with size 'l'
    plot(n,y) %plot the convolution

    %generate value of 'a' as string to set as legend in the graph
    a = num2str(a);
    a = "a="+num2str(a);
    legend(a)

    axis([-12,12,-6,16]) % set axis
    xlabel('n'), ylabel('y[n]') %set labels
    grid on %set grid
end

% the function to return the value of a for the impusle response
function a = get_a_value(ii)
    a_values = [10,3,0,-5]; %all given values of 'a'
    a = a_values(ii); %return value of 'a' to the impulse response
end

% the function to generete the convolution vector y[n] and size of y[n]
function [cnv,L] = linearconvolve(a,b)
    L = length(a)+length(b)-1;
    cnv = zeros(1,L);
    a1=[a,zeros(1,L-length(a))]; % define a new vector of a
    b1=[b,zeros(1,L-length(b))];
    for i=1:L
        c = 0;
        for j=1:i
            c = c + a1(j)*b1(i-j+1);
        end
          cnv(i) = c;
    end
end
