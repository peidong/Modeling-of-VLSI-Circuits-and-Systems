crosstalk
Rd 2 1 300
R1 2 3 100
R2 3 4 150
R3 4 5 250
R4 5 6 300
C1 2 0 2e-14
C2 3 0 5e-14
C3 4 0 1.5e-14
c4 5 0 3.5e-14
C5 6 0 2.3e-14
Cx 7 4 5e-14

VDD 1 0 DC 1V
Vagg 7 0 exp (0 1 0ps 200ps 2000ps 0ps)

.op
.TRAN 1p 1000p

.print all
.plot all
.END
