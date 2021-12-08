FUNCTION spline_evaluate, x, x1, coeff

nx1 = N_ELEMENTS(x1)
y1 = x1 * 0.

FOR i = 0, nx1-1 DO begin
  id = VALUE_LOCATE(x, x1[i])
  y1[i] = coeff.d[id] + $
          coeff.c[id] * (x1[i]-x[id]) + $
          coeff.b[id] * (x1[i]-x[id])^2 + $
          coeff.a[id] * (x1[i]-x[id])^3
ENDFOR

 
          
RETURN, y1
END