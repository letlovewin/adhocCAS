import sys
import math

from inspect import signature

# This is basic calculus stuff.

epsilon = 1/(2**31) # extremely low value

def derivative(f,x):
    return (f(x+epsilon)-f(x))/epsilon

def nthDerivative(f,n,x):
    if n == 1:
        return derivative(f,x)
    return (nthDerivative(f,n-1,x+epsilon)+nthDerivative(f,n-1,x-epsilon))/2*epsilon

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

def intersection(f,g,x_0,max_iter=2**8):
    return newtonRaphson(lambda x: abs(f(x)-g(x)),x_0,max_iter)

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

def shellRevolution(f,a,b):
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

def nthPartialDerivative(f,n,p,wrt): #WARNING!!! THIS ALGORITHM IS QUITE SLOW. IT'S O(2^N). I WROTE THIS NEAR MIDNIGHT AND DID NOT FEEL LIKE DOING ADDITIONAL RESEARCH.
    if n <= 0:
        print("n must be > 0")
        return None
    if wrt <= 0:
        print("wrt must be > 0")
        return None
    if n == 1:
        return nthPartialDerivative(f,n-1,p,wrt)
    newPoint1 = p.copy()
    newPoint1[wrt-1] += epsilon
    newPoint2 = p.copy()
    newPoint2[wrt-1] -= epsilon
    return (nthPartialDerivative(f,n-1,newPoint1)+nthDerivative(f,n-1,x-newPoint2))/2*epsilon

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

def doubleIntegral(f,U):
    # U is a tuple of tuples.
    # If we wanted to integrate over a rectangular region spanning from [0,2] on the x axis and then from [0,4] on the y axis, U = ((0,2),(0,4))

    # OK... again, I wrote a lot of this late at night. I am very sleep deprived and a lot of the papers I'm reading are confusing to me as a result. I won't do anything complicated. This algorithm is an extension of Riemann sums on R to R^n.

    # It uses rectangles of dimensions (b-a)/n * (d-c)/m
    a = U[0][0]
    b = U[0][1]

    c = U[1][0]
    d = U[1][1]

    n = 2**10
    m = 2**10

    delta_x = (b-a)/n
    delta_y = (d-c)/m

    I = 0

    for i in range(a,n+1):
        for j in range(c,m+1):
            
            p1 = f(((a + i * delta_x + a + (i+1) * delta_x)/2,(c + j * delta_y + c + (j+1) * delta_y)/2))
            I += p1*(delta_x)*(delta_y)

    return I
