%Parameter Specification
Rd=300;
Rs=250;
Re=550; 
C1=5e-14;
C2=5.41e-14;
CL=3.89e-14;
Cx=5e-14;
K1=Cx*Rd*Rs*C1;
K2=Rd*Rs*C1*CL*Re*(Cx+C2);
to=200e-12;
%tr = 500e-12;

syms s;
a2=K1/K2;
a1=(Rd+Rs)*Cx/K2;
b2=((C2+Cx)*(Re*CL*(Rd+Rs)+Rd*Rs*C1)+Rd*Re*C1*CL+CL*Rd*Rs*C1)/K2;
b1=((Rd+Rs)*(Cx+C2+CL)+(Re*CL+Rd*C1))/K2;
b0=1/K2;
% Calculate the voltages at each node and the transfer function.
%Vagg = (1-exp(-s*tr))/(s^2*tr);
Vagg = 1/s - 1/(s + 1/to);
display(Vagg);
% V3 = Z3*Vagg/(Z3+1/(s*Cx));
% V4 = V3*(Rr-R3)/Rr;
Vout = (a2*s*s+a1*s)/(s*s*s+b2*s*s+b1*s+b0)*Vagg;
% Hs = Vout/Vagg;

% Hs = simplify(Hs);
% display(Hs);
Vout = simplify(Vout);
digits(1000);
Vout = vpa(Vout);
display(Vout);

%% Plot the output
vout=ilaplace(Vout);
display(vout);
ezplot(1+vout,[0,1000e-12]);
title('Noise Voltage Waveform Considering Crosstalk')
xlabel('t (s)')
ylabel('V (volt)')
%% END