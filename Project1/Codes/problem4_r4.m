mc0 = 1;
mc1 = -7e-11;
m1 = 6e-9;
m2 = -4.36e-19;
m3 = 2.8872e-29;
m4 = -1.88552e-39;

a1 = (m2*m3 - m1*m4) / (m1*m3 - m2*m2);
a2 = (m2*m4 - m3*m3) / (m1*m3 - m2*m2);
p1 = (-a1 + sqrt(a1*a1 - 4*a2))/(2*a2);
p2 = (-a1 - sqrt(a1*a1 - 4*a2))/(2*a2);
k2 = p2 * (mc0/p1 - mc1) / (1/p2 - 1/p1);
k1 = (-mc0 - k2/p2) * p1;

syms s t;
y=ilaplace ((k1/(s-p1) + k2/(s-p2))/s, s, t);
t=0:0.01e-10:1e-10;
% y = (4647675873662383*exp(-(7046493530475741*t)./32768))./56371948243805928 - (2909319386021367*exp(-(8063176603863135*t)./524288))./2687725534621045 + 151512324731213179422899432613341/151512324731213166623240004554760;
y = (38409346964069*exp(-(6256615905320225*t)./65536))./431490752091050 - (8765835262427889*exp(-(2012330311111291*t)./131072))./8049321244445164 + 1736603838794055253056419240567/1736603838794055168215030091100
plot(t,y);
