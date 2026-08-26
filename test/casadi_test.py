import casadi as ca
x = ca.MX.sym('x')
f = (x-3)**2
solver = ca.nlpsol('s','ipopt',{'x':x,'f':f})
sol = solver(x0=0)
print('CasADi+Ipopt OK, x* =', sol['x'])