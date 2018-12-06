# importing the constans , numpy and scipy lib for calculations 
import tables_B1
import tables_B2
import tables_B3
import numpy as np
form sicpy.optimize import fsolve
# This programme is for compression factor calculations at "EN ISO 12213-2:2005" Method.

print('Please Enter the absolute Temperature [K], absolute pressure P [MPa] and the mole fraction of each component "x" as a list')

# Step I : determing B and C*n (n = 13 to 58)

def compCal(T, P, x):
    x = np.array(x)
    N = len(x)
    return x , N

temp1 = [0] * (N-1); temp2 = temp1[:]; B = 0

# Determing of B 

for n in range(18):
    for i in range(N):
        for j in range(N):
            Eij = E_star[i, j] * (E[i] * E[j])**.5
            Gij = G_star[i, j] * (G[i] + G[j])/2
            Bijn_Star[j] = (Gij + 1 - g[n])**g[n] * (Q[i]*Q[j] + 1 - q[n])**q[n] * ((F[i]*F[j])**.5 + 1 - f[n])**f[n] * (S[i]*S[j] + 1 - s[n])**s[n] * (W[i]*W[j] + 1 - w[n])**w[n]
            temp1[j] = x[j] * Bijn_star * E[i] ** (u[n]) * (K[j]) ** (3/4)
        temp2[i] = x[i] * K[i] ** (3/2) * sum(temp1) 
    B += a[n] * T ** (-u[n]) * sum(temp2)

# Determing of Cn_star (n = 13 to 58)

U_cal_term1 = (sum(x * E ** (2.5) )) ** 2
temp0 = 0 ; U_cal_term2 = 0
for i in range(N-1):
    for j in range(i+1 , N):
        temp0 += x[j] * (U[i,j] ** 5 - 1) * (E[j]) ** 2.5
        temp1[i] = temp0
    U_cal_term2 +=  (x[i]) * (E[i]) ** (2.5) * temp1[i])

U_cal_tot = (U_cal_term1 + 2 * U_cal_term2) ** (.2)   # Note that index of each number is one less than the real number

G_cal_term1 = sum(x * G)
for i in range(N-1):
    for j in range(i+1 , N):
        temp0 += x[j] * (G_star[i,j] - 1) * (G[j] + G[i])
        temp1[i] = temp0
    G_cal_term2 +=  (x[i]) * temp1[i])

G_cal_tot = G_cal_term1 + 2 * G_cal_term2

Q_cal = sum(x * Q)

F_cal = sum(x ** 2 * F)

Cn_star = a * (G_cal_tot + 1 -g) ** g * (Q_cal ** 2 + 1 - q) ** q * (F_cal + 1 - f) ** f * U_cal_tot ** u * T ** -u 

# Step III : Calculation of mixture size parameter K and molar density Rho_m

K_cal_term1 = (sum(x * K ** (2.5) )) ** 2
temp0 = 0 ; K_cal_term2 = 0
for i in range(N-1):
    for j in range(i+1 , N):
        temp0 += x[j] * (Kij[i,j] ** 5 - 1) * (K[j]) ** 2.5
        temp1[i] = temp0
    K_cal_term2 +=  (x[i]) * (K[i]) ** (2.5) * temp1[i])
    
K_cal_tot = (K_cal_term1 + 2 * K_cal_term2) ** (.2)   # Note that index of each number is one less than the real number

Rho_m = fsolve(Rho_m * R * T * (1 + B * Rho_m - K**3 * Rho_m * sum(Cn_star[:17]) + sum(Cn_star*(b - c*k*(K**3*Rho_m)**k)*(K**3*Rho_m)**b*exp(-c*(K**3*Rho_m)**k)))) = P , P/(R*T))

# Step IV : The Last is Determing the Z - factor
Z = P/(Rho_m * R * T)
Mw = round(sum(x * M))
Rho = Mw * Rho_m