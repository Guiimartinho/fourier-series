%************************************************************************
% fourier-series-calculate.m
%************************************************************************
% 
% A simple MATLAB script that computes the Fourier series for some kinds of
% signals. Made originally in 2014.
% 
% It shows some basic functions of MATLAB, so may be interesting if you are
% learning how to use it, especially if you are an Electrical Engineering
% student.
% 
% Also a first test with GitHub.
% 
% This script calculates the Fourier Series aproximation of a chosen
% signal: 
%   - Rectangular pulse, with option to choose the ON/OFF ratio 
%   - Exponetial function, with option to choose the value for the decaying
%     constant
%   - Sawtooth function 
%   - Triangular function
% 
% After choosing the function type, you choose how many repetitions.
% 
% Then you choose the form of Fourier Series to be used: 
% - Trigonometric form
% - Compact Trigonometric form
% - Exponential form
% 
% Finally, for the reconstruction, you choose how many coefficients to be 
% used.
% 
% The script then prints the coefficients and plots the function, spectrum
% and phase. The final plot is a comparison of the original and the
% reconstructed signal.
%
%     *****************************************************************
%
%     Written by: Matheus Inoue, State University of New York at New Paltz,
%     2014
%
%     Revised by: Matheus Inoue, Federal University of ABC,  2016
%
%     *****************************************************************


clear all
close all
clc


choice = input ('Select The signal form.\n 1 for Rectangular pulse \n 2 for Exponential \n 3 for Saw-tooth \n 4 for Triagular\n:');

if choice > 4 || choice < 1
    disp('error');
    return;
end

T = input ('Select the period, in seconds: ');
coefs=1:20;%vector used for plotting the coefficients of Fourier Series

%Select the desired signal and the period


n=10000;
dT= T/n;
t=0:dT:T;%vector of time

len=length(t);
tt=t/2;
w0=(2*pi)/T; %omega of function, used in Fourier Series

switch choice
    case 1
        %Rectangular Function
        ratio = input ('Select the ON/OFF ratio: ');
        f= [ones(1,floor(len*ratio)) zeros(1,len - floor(len*ratio))];
        
    case 2
        %Exponential function
        a = input ('Select the value for the coefficient "a": ');
        f=exp(-a.*t);
    case 3
        %Sawtooth function
        f=t;
    case 4
        %Triangular function
        for i=1:length(t)
            if  i< length(t)/2
                f(i)= t(i);
            else
                f(i)= T -t(i);
            end
        end
        
        
end

rep = input ('Choose the number of repetitions: ');
Seriestype = input ('Choose the Fourier Series type. \n 1 for Trigonometric form \n 2 for Compact Trigonometic form \n 3 for Exponential form \n:');

if Seriestype > 3 || choice < 1
    disp('error');
    return;
end

ncoefs = input('Select the number of coefficients to be used: ');

%selection of type of Fourier Series
pow = abs((1/T)*trapz(t,f.^2));%calculate the power of the signal

x=repmat(f,1,rep);
time=linspace(0,rep*T,length(x));%Repetition of selected signal

%plot of the selected signal
plot(time,x,'r')
title('Selected Signal')
xlabel('Time')
ylabel('Amplitude')
grid on;
hold on;
axis([-1,max(time)+1,min(x)-1,max(x)+1])


switch Seriestype
    case 1
        
        a0= (1/T)*trapz(t,f); %DC value
        sumA =0;
        sumB=0;
        sumC=0;
        sumE=0;
        a=0;
        b=0;
        coefs = 0:ncoefs;
        
        if (ncoefs > 0)
        for i=1:ncoefs;
            %Trigonometric Form
            a(i)=(2/T)*trapz(t, f.*cos(i.*w0.*t));
            sumA=a(i)*cos(i.*w0.*t)+sumA;
            b(i)=(2/T)*trapz(t, f.*sin(i.*w0.*t));
            sumB=b(i)*sin(i.*w0.*t)+sumB;
            
        end
        
        a = [a0 a];
        b = [0 b];
                     
        else
        a= a0;
        sumA = zeros(1,length(t));
        b = 0;
        sumB = zeros(1,length(t));
        end
        
        for j=0:ncoefs
            disp(['a',num2str(j),'= ',num2str(a(j+1)),', b',num2str(j),'= ',num2str(b(j+1))])
        end
        
        %plotting the coefficients
        figure;
        stem(coefs,a)
        title('Cosine coefficients')
        xlabel('N')
        ylabel('A_n')
        
        figure;
        stem(coefs,b,'r')
        title('Sine coefficients')
        xlabel('N')
        ylabel('B_n')
        
        
        g= a0 + sumA + sumB;%calculating the signal using fourier series
        form = g;
        
        grep=repmat(g,1,rep);%repeating the signal
        rep = grep;
        
        
    case 2
        
        a0= (1/T)*trapz(t,f);
        sumA =0;
        sumB=0;
        sumC=0;
        sumE=0;
        a=0;
        b=0;
        Cn=0;
        Co=a0;
        theta= 0;
        
        if (ncoefs > 0)
            for i=1:ncoefs;
                %Trigonometric Form
                a(i)=(2/T)*trapz(t, f.*cos(i.*w0.*t));
                sumA=a(i)*cos(i.*w0.*t)+sumA;
                b(i)=(2/T)*trapz(t, f.*sin(i.*w0.*t));
                sumB=b(i)*sin(i.*w0.*t)+sumB;
                
                %Compact Trigonometric Form
                Cn(i)=sqrt(a(i)^2+b(i)^2);
                
                theta(i) = atan2(-b(i),a(i));
                sumC = Cn(i).*cos(i.*w0.*t +theta(i))+sumC;
                
                
            end
                            
            Cn = [a0 Cn];
            theta = [0 (theta*180)/pi];
            
        else
            sumC = zeros(1,length(t));
            Cn(1) = Co;
        end
        
        
        %Displaying the first 20 coefficients and the respective thetas
        coefs = 0:ncoefs;
        for j=0:ncoefs
            disp(['Cn',num2str(j),'= ',num2str(Cn(j+1)),', Theta',num2str(j),' = ',num2str(theta(j+1))])
        end
        
        %Plotting the spectrum and Phase         
        figure;
        grid on;
        stem(coefs,Cn)
        title('Spectrum')
        xlabel('N')
        ylabel('Cn')
        figure;
        stem(coefs,theta, '*r')
        title('Phase')
        xlabel('N')
        ylabel('\theta (degrees)')
        
        
        compform = Co + sumC;%generating signal using fourier series
        form = compform;
        
        crep=repmat(compform,1,rep);
        rep = crep;
        
        
        
    case 3
        
        sumE=0;
        Dn=0;
        j=0+1i;
        c=ncoefs;
        for i=floor(-c):ceil(c)
            Dn(i+ceil(c)+1)=(1/T).*trapz(t,f.*exp(-1*j.*i*w0.*t));
            sumE = Dn(i+ceil(c)+1)*exp(j.*i*w0*t)+ sumE;
        end
        
        %Displaying the coefficients
        for j=-ncoefs:ncoefs
            disp(['D',num2str(j),'= ',num2str(Dn(j+ncoefs+1))])
        end
        
        %Plot of the absolute value and angle
        R = abs(Dn);
        thetaI = angle(Dn);
        thetaI = (thetaI*180)/pi;
        figure;
        d = -ncoefs:ncoefs;
        stem(d,R)
        title('|D_n|')
        xlabel('N')
        figure;
        stem(d,thetaI, '*r')
        title('\angle D_n')
        xlabel('N')
        ylabel('\theta (degrees)')
                
        
        expform = sumE;%generating the signal using Fourier Series
        form = expform;
        
        erep=repmat(expform,1,rep);
        rep = erep;
        
end

minGlob = min(min(real(min(x,rep)),min(abs(rep-x))));
maxGlob = abs(max(max(real(max(x,rep)),max(abs(rep-x)))));

%plotting the signal, the approximation with the Fourier Series and the error
figure;
plot(time,x,'r',time,rep,time,abs(rep-x),'g')
axis([-.5,max(time)+.5,min(minGlob - minGlob*.1, minGlob + minGlob*.1),max(maxGlob - maxGlob*.1, maxGlob + maxGlob*.1)])
legend ('Selected Signal','Reconstructed Signal','Error Signal')
grid on;

%calculating the error
perr=(1/T)*trapz(t,(abs(form-f)).^2);
percentage = (100*(abs(pow-perr)))/pow


switch choice
    case 1
        sigtype = 'Pulse';
    case 2
        sigtype = 'Exponential';
    case 3
        sigtype = 'Saw-tooth';
    case 4
        sigtype = 'Triangular';
    otherwise
        sigtype = 'Unknown';
end
xlabel('Time')
ylabel('Amplitude')
title([sigtype,' signal reconstructed with ',num2str(percentage),'% of the power of the original signal']);