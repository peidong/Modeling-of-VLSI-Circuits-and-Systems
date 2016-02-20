%Parameter Specification
Rd=300;
R1=100;
R2=150;
R3=250;
R4=300; 
C1=2e-14;
C2=5e-14;
C3=1.5e-14;
C4=3.5e-14;
C5=2.3e-14; 
Cx=5e-14;
tr = 500e-12;
% Calculate the Impedances in s-domain.
syms s;
Z1 = 1/(1/Rd+s*C1);
Z2 = 1/(1/(Z1+R1)+s*C2);
Rr = 1/(1/(R4+1/(s*C5))+s*C4)+R3;
Z3 = 1/(1/(Z2+R2)+s*C3+1/Rr);
% Calculate the voltages at each node and the transfer function.
Vagg = (1-exp(-s*tr))/(s^2*tr);
display(Vagg);
V3 = Z3*Vagg/(Z3+1/(s*Cx));
V4 = V3*(Rr-R3)/Rr;
Vout = V4*1/(s*C5)/(R4+1/(s*C5));
Hs = Vout/Vagg;

Hs = simplify(Hs);
display(Hs);
Vout = simplify(Vout);
digits(10);
Vout = vpa(Vout);
display(Vout);

%% Plot the output
vout=ilaplace(Vout);
display(vout);
ezplot(vout,[0,800e-12]);
title('Noise Voltage Waveform Considering Crosstalk')
xlabel('t (s)')
ylabel('V (volt)')
%% END