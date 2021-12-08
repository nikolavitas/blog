FUNCTION csc_interpolation, x, y, yprim = yprim

; Define the variables
n = SIZE(x, /N_ELEMENTS)-1
h = DBLARR(n+1)
w = DBLARR(n+1)
q = DBLARR(n+1)
d = DBLARR(n+1)
c = DBLARR(n+1)
b = DBLARR(n+1)
a = DBLARR(n+1)
ctemp = DBLARR(n+1)

; Compute x-differences
h[0] = x[1] - x[0]
h[n-1] = x[n] - x[n-1]

; Compute the first derivative using the IDL intrinsic function.
IF ~KEYWORD_SET(yprim) THEN yprim = DERIV(x, y)

w[0] = 2*h[0]
w[n] = 2*h[n-1]

q[0] = (3.D/h[0])*(y[1] - y[0]) - 3*yprim[0]
q[n] = 3*yprim[n] - (3.D/h[n-1])*(y[n] - y[n-1])


; Compute H, W and Q
FOR i = 1, n-1 DO BEGIN
  h[i] = x[i+1] - x[i]
  w[i] = 2*(x[i+1] - x[i-1])
  q[i] = 3*(y[i+1] - y[i])/h[i] - 3*(y[i] - y[i-1])/h[i-1]
ENDFOR

; Gaussian elimination to get W' and Q'
FOR i = 1, n DO BEGIN
  w[i] = w[i] - h[i-1]^2/w[i-1]
  q[i] = q[i] - q[i-1]*h[i-1]/w[i-1]
ENDFOR

; Solve the system by succesive backward subsitutes

;ctemp[n] = 0
ctemp[n] = q[n]/w[n]

FOR i = 1, n DO ctemp[n-i] = (q[n-i] - h[n-i]*ctemp[n-i+1])/w[n-i]

; 0-boundary for the spline parameters
d[0] = (ctemp[1]-ctemp[0])/(3*h[0])
c[0] = ctemp[0]
b[0] = yprim[0]
b[n] = yprim[n]
a[0] = y[0]

; Compute all the parameters.
FOR i = 1, n-1 DO BEGIN
  d[i] = (ctemp[i+1] - ctemp[i])/(3*h[i])
  c[i] = ctemp[i]
  b[i] = (ctemp[i] + ctemp[i-1])*h[i-1] + b[i-1]
  a[i] = y[i]
ENDFOR

result = {a:a, b:b, c:c, d:d}

RETURN, result
END