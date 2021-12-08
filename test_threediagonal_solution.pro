FUNCTION test_threediagonal_solution, a, b, c, d, x

n = N_ELEMENTS(d)

d0 = FLTARR(n)

d0[0] = b[0]*x[0] + c[0]*x[1]
FOR i = 1, n-2 DO d0[i] = a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1]
d0[n-1] = a[n-1]*x[n-2] + b[n-1]*x[n-1] 

RETURN, d0
END