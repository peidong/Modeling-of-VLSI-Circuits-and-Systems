*struct
V1 6 0 pulse(0 1 0 0.001ps)
R1 6 1 2m
R2 1 2 2m
R3 2 3 3m
R4 2 4 3m
R5 4 5 4m
C1 1 0 2n
C2 2 0 2n
C3 3 0 4n
C4 4 0 4n
C5 5 0 2n
.op
.tran 0.1ps 500ps
.print all
.plot all
.end
