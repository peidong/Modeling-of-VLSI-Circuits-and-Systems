
VDD s 0 PULSE(0 1 0 0.001ps)
R1 s 1 2m
R2 1 2 2m
R3 2 3 3m
R4 2 n 3m
R5 n 5 4m
C1 1 0 2n
C2 2 0 2n
C3 3 0 4n
C4 n 0 4n
C5 5 0 2n

.op
.TRAN 0.1ps 500ps
.print all
.plot all
.END
