R = QQ[s,t]
S = R[y12,y13,y14,y23,y24,y34]
z1 = y12 + y13 + y14
z2 = -y12 + y23 + y24
z3 = -y13 - y23 + y34
z4 = -y14 - y24 - y34
I = ideal(y12*y13 - y12*y23 + y13*y23,
         y12*y14 - y12*y24 + y14*y24,
         y13*y14 - y13*y34 + y14*y34,
         y23*y24 - y23*y34 + y24*y34,
         -t * z1*z2 + (s-t)*y12^2,
         -t * z1*z3 + (s-t)*y13^2,
         -t * z1*z4 + (s-t)*y14^2,
         -t * z2*z3 + (s-t)*y23^2,
         -t * z2*z4 + (s-t)*y24^2,
         -t * z3*z4 + (s-t)*y34^2)
T = S/(I + truncate(4, ideal vars S))
M = image basis(3,T)
N = pushForward(map(T,R), M)
A = annihilator N
primaryDecomposition A
