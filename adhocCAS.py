import sys
import math

from inspect import signature

# This is basic calculus stuff.

epsilon = 1/(2**31) # extremely low value

def derivative(f,x):
    return (f(x+epsilon)-f(x))/epsilon

def newtonRaphson(f,x_0,max_iter): # Find a root of f given an initial guess x_0.
    xn = x_0
    for k in range(max_iter + 1):
        fxn = f(xn)
        if(abs(fxn)<=epsilon):
            return xn

        Dx = derivative(f,xn)
        if Dx == 0:
            return None

        xn = xn - fxn/Dx
    return None

def riemannSum(f,a,b): # Find the area under the curve of f from a to b.
    if a == b:
        return 0
    if a > b:
        return -1 * riemannSum(b,a)
    S = 0
    n = 2**16
    for k in range(0,n):
        S += f(a+k*((b-a)/n))
    return S*(b-a)/n

def arcLength(f,a,b):
    return riemannSum(lambda x: math.sqrt(1+derivative(f,x)**2),a,b)

def areaBounded(f,g,a,b):
    return riemannSum(lambda x: abs(f(x)-g(x)),a,b)

def volumeOfRevolution(f,g,a,b): # Find the volume of a curve f revolved around an axis g from a to b.
    return riemannSum(lambda x:abs((f(x)-g(x)))**2,a,b)*3.14159

def shellRevolution(f,a,b): # Find the volume of a curve f revolved around an axis g from a to b. This is a shell
    return 3.14159*riemannSum(lambda x:x*f(x),a,b)

# Multivariate functions.

def partialDerivative(f,p,wrt): # The function you pass here should take a list/tuple, e.g. lambda x:x[0]+x[1] for some function f(x_1,x_2) = x_1 + x_2
    # A sample call of this function is partialDerivative(lambda x: x[0]**2 + x[1]**2,[2,2],2), giving us the partial derivative of our function wrt y at the point (2,2).
    if wrt <= 0:
        print("wrt must be > 0")
        return None
    newPoint = p.copy()
    newPoint[wrt-1] += epsilon
    return (f(newPoint) - f(p))/epsilon

def gradient(f,p):
    r = []
    for i in range(len(p)):
        r.append(partialDerivative(f,p,i+1))
    return r

def parametricSub(f,g,t):
    S = 0
    k = []
    for i in range(len(g)):
        k.append(g[i](t))
    return f(k)

def Ds(c,x):
    S = 0
    for i in range(len(c)):
        S += derivative(c[i],x)**2
    return math.sqrt(S)

def lineIntegral(f,c,a,b):
    # f is the surface that you'll be integrating over.
    # c is the parametric line that you'll be integrating under, let c be a tuple of functions, e.g. c = [lambda x: math.sin(x),lambda x: math.cos(x)]
    # Say we wanted to integrate over xy^4 with the parameterization (4cos(t),4sin(t)), from -pi/2 to pi/2. Then we'd call lineIntegral(lambda x: x[0]*x[1]**4, (lambda t: 4*math.cos(t),lambda t: 4*math.sin(t)), -1.57079,1.57079)
    return riemannSum(lambda x: parametricSub(f,c,x)*Ds(c,x),a,b)