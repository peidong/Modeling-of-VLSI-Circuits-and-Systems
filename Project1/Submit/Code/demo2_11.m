clear;

%load G,C matrices
a=load('GC8.txt');
% a=load('GC9.txt');
G=sparse(a(:,1),a(:,2),a(:,3));
C=sparse(a(:,1),a(:,2),a(:,4));
kv=3;   %number of voltage sourse

%define matched moments, size of original system, number of input sources
q=5;         %define matched moments
N=size(G,1);
gmin=0; %perturbation to make reduced G nonsigular

%define start frequency, end frequency, and step
fspan=[1e3 1e5 1e7 1e9];      %span vectors at the frequency of 1MHz
fe=1e9;

%set parameters for time domain response calculation
h=2e-9;       %step for back-ward time domain simulation
inputno=8000;     %the number of time steps that have input
outputno=inputno;  %the number of time steps for output we consider

%Generate L
L=[zeros(N,1)]; 
L(8)=1;

%Generate B
B=zeros(size(G,1),kv);
B(6,1)=-1;
B(10,2)=-1;
B(11,3)=-1;

%Generate U                    
%
%Voltage source 1      0V    --------------------   
%Voltage source 2       0.5V*sin(2*pi*1MHz*t)       
%Voltage source 3      -0.5V*sin(2*pi*1MHz*t)
U=[0*[h:h:h*inputno];0.5*sin(2*pi*1e6*[h:h:h*inputno]);-0.5*sin(2*pi*1e6*[h:h:h*inputno])];
Us=fft(U,inputno*1000,2);
f=1/h/2*linspace(0,1,inputno*1000/2);

G=G+gmin*eye(length(G));

%Prima reduction
fprintf('\n\n\n\n');
fprintf('G,C,B,U,L matrices have been generated.');

fprintf('\n');
fprintf('Prima begins:\n');
tic
[Gr,Cr,Br,Lr,V]=prima(G,C,B,L,q,2*pi*fspan,gmin);
toc
fprintf('Prima done!\n');

%calculate original time domain response
fprintf('\n');
fprintf('Calculate original time domain response:\n');
tic
vo=-1e-9*ones(N,1); 
A=G+1/h*C;
[LA,UA]=lu(A);

for j=1:outputno
    if (j<=inputno)
        b=1/h*C*vo(:,end) + B*U(:,j);
    else
        b=1/h*C*vo(:,end);
    end
    vo=[vo, UA\(LA\b)];
end

To=L'*vo;

toc
fprintf('Original time domain response done!\n');

%calculate reduced time domain response
fprintf('\n');
fprintf('Calculate reduced time domain response:\n');
tic
vr=-1e-9*ones(size(Gr,1),1); 
Ar=Gr+1/h*Cr;
[LAr,UAr]=lu(Ar);

for j=1:outputno
    if (j<=inputno)
        br=1/h*Cr*vr(:,end) + Br*U(:,j);
    else
        br=1/h*Cr*vr(:,end);
    end
    vr=[vr, UAr\(LAr\br)];
end

Tr=Lr'*vr;

toc
fprintf('Reduced time domain response done!\n');
 
%calculate original frequency domain response
fprintf('\n');
fprintf('Calculate original frequency response:\n');
tic
cnt=0;
bb=size(B);
Fo=[];
for row_idx=0:0.02:log10(length(f))
    s=2*pi*f(floor(10^(row_idx)))*1i;
    z=L'*((G+s*C)\(B*Us(:,floor(10^(row_idx)))));
    Fo=[Fo,z];
end
toc
fprintf('Original frequency response done!\n');

%calculate reduced frequency domain response
fprintf('\n');
fprintf('Calculate reduced frequency response:\n');
tic
cnt=0;
bb=size(Br);
Fr=[];
for row_idx=0:0.02:log10(length(f))
    s=2*pi*f(floor(10^(row_idx)))*1i;
    z=Lr'*((Gr+s*Cr)\(Br*Us(:,floor(10^(row_idx)))));
    Fr=[Fr,z];
end
toc
fprintf('Reduced frequency response done!\n');

%plot
figure
subplot(1,3,1);
semilogx(f(floor(10.^[0:0.02:log10(length(f))])),10*log10(abs(Fo)),'r');
hold on;
semilogx(f(floor(10.^[0:0.02:log10(length(f))])),10*log10(abs(Fr)),'b');
legend('Original','Reduced');
title('frequency response, magnitude');
hold off;

subplot(1,3,2);
semilogx(f(floor(10.^[0:0.02:log10(length(f))])),angle(Fo),'r');
hold on;
semilogx(f(floor(10.^[0:0.02:log10(length(f))])),angle(Fr),'b');
legend('Original','Reduced');
title('frequency response, phase');
hold off;

subplot(1,3,3);
plot(h*[1:length(To)],To,'r');
hold on;
plot(h*[1:length(Tr)],Tr,'b');
legend('Original','Reduced');
title('time response');
hold off;

% Calculate inpulse response
fprintf('\nCalculate original impulse response:\n');
tic
cnt=0;
bb=size(B);
Fo2=L'*(G\(B*ones(bb(2),1)));
for row_idx=0:0.02:log10(fe)
    s=2*pi*(10^row_idx)*1i;
    z=L'*((G+s*C)\(B*[0;0.5;-0.5]));
    Fo2=[Fo2,z];
end
toc
fprintf('Original impulse response done!\n');

fprintf('\nCalculate reduced impulse response:\n');
tic
cnt=0;
bb=size(Br);
Fr2=Lr'*((Gr\(Br*ones(bb(2),1))));
for row_idx=0:0.02:log10(fe)
    s=2*pi*(10^row_idx)*1i;
    z=Lr'*((Gr+s*Cr)\(Br*[0;0.5;-0.5]));
    Fr2=[Fr2,z];
end
toc
fprintf('Reduced impulse response done!\n\n');

figure;
subplot(1,2,1);
semilogx([0,10.^[0:0.02:log10(fe)]],abs(Fo2),'r');
hold on;
semilogx([0,10.^[0:0.02:log10(fe)]],abs(Fr2),'b');
legend('Original','Reduced');
title('frequency response, magnitude');
hold off;

subplot(1,2,2);
semilogx([0,10.^[0:0.02:log10(fe)]],angle(Fo2),'r');
hold on;
semilogx([0,10.^[0:0.02:log10(fe)]],angle(Fr2),'b');
legend('Original','Reduced');
title('frequency response, phase');
hold off;
