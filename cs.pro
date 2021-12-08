FUNCTION cs, x, y, x0, spline = spline, output = output, dy = dy, d2y = d2y, d3y = d3y

nx = N_ELEMENTS(x)
nx0 = N_ELEMENTS(x0)

IF ~KEYWORD_SET(spline) THEN spline = 'natural'
IF ~KEYWORD_SET(output) THEN output = 'standard'

IF spline EQ 'natural' THEN coeff = CSN_INTERPOLATION(x, y)
IF spline EQ 'clamped' THEN coeff = CSC_INTERPOLATION(x, y)

a = coeff.a
b = coeff.b
c = coeff.c
d = coeff.d

IF output EQ 'standard' THEN BEGIN
  y0 = DBLARR(nx0)
  FOR i = 0, N_ELEMENTS(x0)-1 DO BEGIN
    j = VALUE_LOCATE(x, x0[i]) 
    IF (j GE 0) AND (j LT nx-1) THEN BEGIN
      y0[i] = a[j] + b[j]*(x0[i]-x[j]) + c[j]*(x0[i]-x[j])^2 + d[j]*(x0[i]-x[j])^3
      dy[i] = b[j] + 2*c[j]*(x0[i]-x[j]) + 3*d[j]*(x0[i]-x[j])^2
      d2y[i] = 2*c[j] + 6*d[j]*(x0[i]-x[j])
      d3y[i] = 6*d[j]
    ENDIF
    IF (j GE 0) AND (j EQ nx-1) AND (x[j] EQ x0[i]) THEN BEGIN
      y0[i] = y[j]
      dy[i] = b[j-1] + 2*c[j-1]*(x0[i]-x[j-1]) + 3*d[j-1]*(x0[i]-x[j-1])^2
      d2y[i] = 2*c[j-1] + 6*d[j-1]*(x0[i]-x[j-1])
      d3y[i] = 6*d[j-1]
    ENDIF
  ENDFOR
ENDIF  

IF output EQ 'full' THEN BEGIN
  y0 = DBLARR(nx, nx0)+missing
  FOR j = 0, nx-2 DO BEGIN
    FOR i = 0, nx0-1 DO BEGIN
      y0[j, i] = coeff.a[j] + coeff.b[j]*(x0[i] - x[j]) + $
                 coeff.c[j]*(x0[i] - x[j])^2 + coeff.d[j]*(x0[i] - x[j])^2
    ENDFOR
  ENDFOR
ENDIF  


RETURN, y0
END