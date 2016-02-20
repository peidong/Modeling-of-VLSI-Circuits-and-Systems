ms0 = 1;
ms1 = -7e-11;
m1 = 1.4000e-8;
m2 = -8.5200e-19;
m3 = 5.5128e-29;
m4 = -3.5833e-39;
 
a1 = (m2*m3 - m1*m4) / (m1*m3 - m2*m2);
a2 = (m2*m4 - m3*m3) / (m1*m3 - m2*m2);
p1 = (-a1 + sqrt(a1*a1 - 4*a2))/(2*a2);
p2 = (-a1 - sqrt(a1*a1 - 4*a2))/(2*a2);
k2 = p2 * (ms0/p1 - ms1) / (1/p2 - 1/p1);
k1 = (-ms0 - k2/p2) * p1;
 
syms s t;
y=ilaplace ((k1/(s-p1) + k2/(s-p2))/s, s, t);
t=0:0.01e-10:1e-10; 
y = (4647675873662383*exp(-(7046493530475741*t)./32768))./56371948243805928 - (2909319386021367*exp(-(8063176603863135*t)./524288))./2687725534621045 + 151512324731213179422899432613341/151512324731213166623240004554760;
plot(t,y);
