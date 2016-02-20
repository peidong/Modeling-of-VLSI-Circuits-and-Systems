% Copyright (c) University of California, Los Angeles. All rights reserved.
%                EE201C: Modeling of VLSI Circuits and Systems
%                             2016 Winter         (due on 02/01/2016)
%
%

%% - Step 1: discretization -

% - three wires -
% - Unit: um (1e-6 m) -
for i=1:3
    wire{i}.width=9e-6;
    wire{i}.thickness=6e-6;
    wire{i}.length= 9000e-6;
    wire{i}.x=(i-1)*24*1e-6;
    wire{i}.y=0;
end

% - discretize each wire into two filaments -
for i=1:6
    filament{i}.width=4.5e-6;
    filament{i}.thickness=6e-6;
    filament{i}.length=9000e-6;
    filament{i}.y=0;
    if mod(i,2)~=0
        filament{i}.x=(i-1)/2*24*1e-6;
    else
        filament{i}.x=((i-2)/2*24+4.5)*1e-6;
    end
end

%% - Step 2: Inductance Calculation [without ground] -
ind_filament = zeros(6,6);
% - (2.1) self-Inductance of filaments -
u = 1.256e-6;
for i = 1:6
    % TASK: use equation to compute self-inductance for each filament
    ind_filament(i,i) = (u*filament{i}.length/(2*pi))*(log(filament{i}.length*2/(filament{i}.width+filament{i}.thickness))+0.5+(filament{i}.width+filament{i}.thickness)/(4*filament{i}.length));
end

% - (2.2) mutual-Inductance of filaments -
for i=1:6
    for j=1:6
        % TASK: use equation to compute mutual-inductance between each pair of filaments
        if(i~=j)
        ind_filament(i,j) = (u*filament{i}.length/(2*pi))*(log(filament{i}.length*2/(abs(filament{i}.x-filament{j}.x))) - 1 + (abs(filament{i}.x-filament{j}.x))/filament{i}.length);
        end
    end
end

% - (2.3) calculate wire inductance -
ind_wire = zeros(3,3);
for i=1:3
    for j=1:3
        if (i==j)
        %TASK: use equation to calculate self-inductance of each wire
        ind_wire(i,j) = ind_filament(2*i-1,2*j-1)+ind_filament(2*i,2*j)+ind_filament(2*i-1,2*j)+ind_filament(2*i,2*j-1);%rref
        else
        %TASK: use equatoin to calculate mutual-inductance of each pair of wires
        ind_wire(i,j) = ind_filament((i-1)*2+1,(j-1)*2+1)+ind_filament((i-1)*2+1,(j-1)*2+2)+ind_filament((i-1)*2+2,(j-1)*2+1)+ind_filament((i-1)*2+2,(j-1)*2+2);
        end
    end
end
K1 = ind_wire(1,2)/ind_wire(1,1);
K2 = ind_wire(1,3)/ind_wire(1,1);
K3 = ind_wire(2,3)/ind_wire(3,3);


%% - Step 3: Capcitance Calculation [with ground]  - 

% - (3.1) capacitance of signal wire -

% TASK: use equation to calculate the capacitance between signal wire (middle one) and
% ground

w=9e-6;
t=6e-6;
l=9000e-6;
s=15e-6;
h=10e-6;
e=8.85e-12;

% TASK: use equation to calculate the coupling capacitance
C3  = e*(1.15*(w/h)+2.8*(t/h)^0.222)*l; %single wire
C2  = e*(0.03*(w/h)+0.83*(t/h)-0.07*(t/h)^0.222)*(s/h)^-1.34*l; %coupling capacitance

C4  = C2;

C31 = C3/2;
C32 = C3/2;


% - (3.2) coupling capacitance of edge wires-

% TASK: use equation to calculate capcitance between edge wire and ground
%C1 = C5;
C1  = C2+C3;
C5  = C1;
C11 = C1/2;
C12 = C1/2;
C51 = C5/2;
C52 = C5/2;



%% - Step 4: Resistance Calculation
r_wire = zeros(3,1);
for i=1:3
    % TASK: calculate R for each wire: r_wire(i)
    r_wire(i) = 0.0175*l/w/t/10e5;
end
    
%% - Step 5: Generate RC and RCL netlist for SPICE -
fid1 = fopen('rc.sp', 'w');
fprintf(fid1, '*This is RC circuit\n');
%TASK: write all R and C component into the Spice Input file
% fprintf(fid1, 'VDD 1 0 PULSE(0 1 0 10ps)\n');
fprintf(fid1, 'VDD 1 0 PULSE(0 1 0 20ns)\n');
fprintf(fid1, 'R1 3 4 %e\n',r_wire(1));
fprintf(fid1, 'R2 1 2 %e\n',r_wire(2));
fprintf(fid1, 'R3 5 6 %e\n',r_wire(3));
fprintf(fid1, 'C11 3 0 %e\n',C11);
fprintf(fid1, 'C12 4 0 %e\n',C12);
fprintf(fid1, 'C31 1 0 %e\n',C31);
fprintf(fid1, 'C32 2 0 %e\n',C32);
fprintf(fid1, 'C51 5 0 %e\n',C51);
fprintf(fid1, 'C52 6 0 %e\n',C52);
fprintf(fid1, 'C2 4 2 %e\n',C2);
fprintf(fid1, 'C4 6 2 %e\n',C4);
fprintf(fid1, '.op\n');
% fprintf(fid1, '.TRAN 0.1ps 30ps\n');
fprintf(fid1, '.TRAN 0.1ps 50ns\n');
fprintf(fid1, '.print all\n');
fprintf(fid1, '.plot all\n');
fprintf(fid1, '.END\n');
fclose(fid1);


fid2 = fopen('rlc.sp', 'w');
fprintf(fid2, '* This is RLC circuit\n');
% TASK: write all R, L and C component into the Spice Input file
fprintf(fid2, 'VDD 1 0 PULSE(0 1 0 10ps)\n');
%fprintf(fid2, 'VDD 1 0 PULSE(0 1 0 20ps)\n');
fprintf(fid1, 'R1 4 5 %e\n',r_wire(1));
fprintf(fid1, 'R2 1 2 %e\n',r_wire(2));
fprintf(fid1, 'R3 7 8 %e\n',r_wire(3));
fprintf(fid1, 'C11 4 0 %e\n',C11);
fprintf(fid1, 'C12 6 0 %e\n',C12);
fprintf(fid1, 'C31 1 0 %e\n',C31);
fprintf(fid1, 'C32 3 0 %e\n',C32);
fprintf(fid1, 'C51 7 0 %e\n',C51);
fprintf(fid1, 'C52 9 0 %e\n',C52);
fprintf(fid1, 'C2 6 3 %e\n',C2);
fprintf(fid1, 'C4 9 3 %e\n',C4);
fprintf(fid1, 'L1 5 6 %e\n',ind_wire(1,1));
fprintf(fid1, 'L2 2 3 %e\n',ind_wire(1,1));
fprintf(fid1, 'L3 8 9 %e\n',ind_wire(1,1));
fprintf(fid1, 'K1 L1 L2 %e\n',K1);
fprintf(fid1, 'K2 L1 L3 %e\n',K2);
fprintf(fid1, 'K3 L2 L3 %e\n',K3);
fprintf(fid2, '.op\n');
fprintf(fid2, '.TRAN 0.1ps 5ns\n');
fprintf(fid2, '.print all\n');
fprintf(fid2, '.plot all\n');
fprintf(fid2, '.END\n');
fclose(fid2);


%% - Run SPICE to compare waveforms from different models -

% - End -
