clear all
close all
clc


choice = input ('Select The signal form.\n 1 for Rectangular \n 2 for Exponential \n 3 for Saw-tooth \n 4 for Triagular\n:');
if choice > 4 || choice < 1
   disp('error');
   return;
end

T = input ('Select the period, in seconds: ');
coefs=1:20;%vector used for plotting the coefficients of Fourier Series

%Select the desired signal and the period


n=1000;
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
         f=t;
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
%selection of type of Fourier Series
pow = (1/T)*trapz(t,f);%calculate the power of the signal

x=repmat(f,1,rep);
time=linspace(0,rep*T,length(x));%Repetition of selected signal

plot(time,x,'r')
title('Selected Signal')
 xlabel('Time')
 ylabel('Amplitude')
 hold on
 axis([-1,max(time)+1,min(x)-1,max(x)+1])
%plot of the selected signal


%Fourier Series Coefficients

a0= (1/T)*trapz(t,f);%DC value
sumA =0;
sumB=0;
sumC=0;
sumE=0;
for i=1:20;
    %Trigonometric Form
   a(i)=(2/T)*trapz(t, f.*cos(i.*w0.*t));% An coefficient
   sumA=a(i)*cos(i.*w0.*t)+sumA; %Summation of An terms
   b(i)=(2/T)*trapz(t, f.*sin(i.*w0.*t)); %Bn coefficient
    sumB=b(i)*sin(i.*w0.*t)+sumB; %Summation of Bn terms
    
    %Compact Trigonometric Form
    Cn(i)=sqrt(a(i)^2+b(i)^2);%Cn coefficient
    theta(i) = atan2(-b(i),a(i));%Phase
    sumC = Cn(i).*cos(i.*w0.*t +theta(i))+sumC;%summation of Cn
    
          
end
 
%For the Exponential form
for i=-10:10
    Dn(i+11)=(1/T).*trapz(t,f.*exp(-1*j.*i*w0.*t));%Dn coefficient
    sumE = Dn(i+11)*exp(j*i*w0*t)+sumE;%summation of Dn
end

if Seriestype ==1
  
 
    for j=1:20
        disp(['a',num2str(j),'= ',num2str(a(j)),', b',num2str(j),'= ',num2str(b(j))])
    end
%Displaying the first 20 coefficients

figure
plot(coefs,Cn,'o')
title('Spectrum')
xlabel('N')
ylabel('Cn')

figure
plot(coefs,theta,'*')
title('Phase')
xlabel('N')
ylabel('Theta')
%plotting Spectrum and Phase


    elseif Seriestype ==2
    %Compact Form
    Co=a0;


    for j=1:20
        disp(['Cn',num2str(j),'= ',num2str(Cn(j)),' Theta',num2str(j),' = ',num2str(theta(j))])
    end
%Displaying the first 20 coefficients and the respective thetas

figure
grid on
plot(coefs,Cn,'o')
title('Spectrum')
xlabel('N')
ylabel('Cn')
figure
plot(coefs,theta, '*')
title('Phase')
xlabel('N')
ylabel('Theta')

%plotting the Spectrum and Phase



elseif Seriestype ==3

%Exponential Form

for j=-10:9
    disp(['D',num2str(j),'= ',num2str(Dn(j+11))])
end
%Displaying the first 20 coefficients

R = abs(Dn);
thetaI = angle(Dn);
figure
d=1:21;
plot(d,R, 'o')
title('Spectrum')
xlabel('N')
ylabel('Dn')
figure
plot(d,thetaI, '*')
   title('Phase')
xlabel('N')
ylabel('Theta')
%plotting Absolute value of Dn and its phase
        end


if Seriestype == 1
    
    ncoefs = input('Select the number of coefficients to be used: ');

a0= (1/T)*trapz(t,f);
sumA =0;
sumB=0;
sumC=0;
sumE=0;
a=0;
b=0;

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





 
 
g= a0+sumA + sumB;%calculating the signal using fourier series
grep=repmat(g,1,rep);%repeating the signal
time=linspace(0,rep*T,length(x));


figure
plot (time,x,'r',time, grep,time,abs(grep-x),'g')
legend ('Selected Signal','Reconstructed Signal','Error Signal')
axis([-1,max(time)+1,min(x)-1,max(x)+1])
%plotting the graph of selected signal, the signal generated and the error
perr=(1/T)*trapz(t,(abs(g-f)).^2);
   percentage = (100*(pow-perr))/pow
%calculus of the percentage

    
elseif Seriestype ==2
    
ncoefs = input('Select the number of coefficients to be used: ');
a0= (1/T)*trapz(t,f);
sumA =0;
sumB=0;
sumC=0;
sumE=0;
a=0;
b=0;
Cn=0;

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
compform = Co + sumC;%generating signal using fourier series


 

crep=repmat(compform,1,rep);

figure
plot (time,x,'r',time,crep,time,abs(crep-x),'g')
axis([-1,max(time)+1,min(x)-1,max(x)+1])
legend ('Selected Signal','Reconstructed Signal','Error Signal')
%plotting the signal, the generated one and the error
perr=(1/T)*trapz(t,(abs(compform-f)).^2);
   percentage = (100*(pow-perr))/pow
%calculating the power of signal





elseif Seriestype ==3

ncoefs = input('Select the number of coefficients to be used: ');
sumE=0;
Dn=0;
j=0+1i;
c=ncoefs/2;
for i=floor(-c):floor(c)
    Dn(i+11)=(1/T).*trapz(t,f.*exp(-1*j.*i*w0.*t));
    sumE = Dn(i+      11)*exp(j*i*w0*t)+sumE;
end

expform = sumE;%generating the signal using Fourier Series


erep=repmat(expform,1,rep);

figure
plot(time,x,'r',time,erep,time,abs(erep-x),'g')
% axis([-1,max(time)+1,min(x)-1,max(x)+1])
legend ('Selected Signal','Reconstructed Signal','Error Signal')
grid on
%plotting the signal, the one with Fourier Series and the error
   perr=(1/T)*trapz(t,(abs(expform-f)).^2);
   percentage = (100*(pow-perr))/pow
%calculating the error

end