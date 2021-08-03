\ENCODING=ISO-8859-1
\Problem name: ILOG.CPLEX

Minimize
 obj: - 1.5 x_0 - 4 x_1 + theta
Subject To
 c1: 2 x_0 + 3 x_1 <= 4
 c2: x_0 >= 1
 c3: x_1 <= 0
 c5: x_0  = 0
Lazy Constraints
 l1: x_0 >= 1
 l2: x_1 <= 0
Bounds
 0 <= x_0 <= 1
 0 <= x_1 <= 1
-99 <= theta <= 99
Binaries
 x_0  x_1 
End
