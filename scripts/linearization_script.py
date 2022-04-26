import sympy

x = sympy.Symbol('x')
y = sympy.Symbol('h')
f = x * x
df = sympy.diff(f,'x')
dff = sympy.diff(f,'x').subs(df,'dx')
print(df)
print(dff)
#ode = f(x).diff(x, 2) - sin(f(x))
#ode_linear = series(ode.subs(f(x), y), y, 0, 2).removeO().subs(y, f(x))
#ode_cubic = series(ode.subs(f(x), y), y, 0, 4).removeO().subs(y, f(x))