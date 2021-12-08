; This file is not meant to be shared. It's an illustration of
; how to run the csn_interpolation routine.

PRO cs_example1

!p.background = GETCOLOR('white')
!p.color = GETCOLOR('black')

; The 'unknown' function, finely sampled at 201 point at the [-1, 1] interval
nxa = 201
xa = DINDGEN(nxa)/100-1.
ya = CS_FUNCTION(xa)

; The 'measurement' sample coarsly the function at 9 points
x = DOUBLE([-1., -0.8, -0.6, -0.45,  0.0, 0.3, 0.5, 0.6, 1.])
; x = RANDOMU(seed, 28)*2. - 1.
; x = [-1, x[SORT(x)], 1]
nx = N_ELEMENTS(x)
y = CS_FUNCTION(x)

; Compute the approximative function S and evaluate it at the coarse grid.
x0 = xa
nx0 = nxa
; y0 = CS(x, y, x0, spline = 'natural', dy = dy, d2y = d2y, d3y = d3y)

; Do the interpolation
coeff = CSN_INTERPOLATION(x, y)

a = coeff.a
b = coeff.b
c = coeff.c
d = coeff.d

y0 = DBLARR(nx0)
dy0 = DBLARR(nx0)
d2y0 = DBLARR(nx0)
FOR i = 0, nx0-1 DO BEGIN
  j = VALUE_LOCATE(x, x0[i]) 
  IF (j GE 0) AND (j LT nx-1) THEN BEGIN
    y0[i] = a[j] + b[j]*(x0[i]-x[j]) + c[j]*(x0[i]-x[j])^2 + d[j]*(x0[i]-x[j])^3
    dy0[i] = b[j] + 2*c[j]*(x0[i]-x[j]) + 3*d[j]*(x0[i]-x[j])^2
    d2y0[i] = 2*c[j] + 6*d[j]*(x0[i]-x[j])
  ENDIF
  IF (j GE 0) AND (j EQ nx-1) THEN BEGIN
    IF (x[j] EQ x0[i]) THEN BEGIN
      y0[i] = y[j]
      dy0[i] = b[j-1] + 2*c[j-1]*(x0[i]-x[j-1]) + 3*d[j-1]*(x0[i]-x[j-1])^2
      d2y0[i] = 2*c[j-1] + 6*d[j-1]*(x0[i]-x[j-1])
    ENDIF
  ENDIF
ENDFOR

name = 'noborder'
; 
; ; Fig.1 The exact "unknown" function is shown in black. The locations of
; ;       the "measurements" are labeled with pluses. The interpolation
; ;       function is shown in red. 
WINDOW, 0, xsi = 500, ysi = 400
PLOT, xa, ya, xr = [-1., 1], /xs, /ys, xtit = 'x', ytit = 'y', tit = 'S', thick = 2, $
      charsi = 1.4, yr = [-0.5, 0.5]
OPLOT, x, y, psym = 1, symsi = 1.4
OPLOT, x0, y0, col = getcolor('red'), thick = 2
WRITE_PNG, 'example_CSN_0'+name+'.png', tvrd(/true)

WINDOW, 0, xsi = 500, ysi = 400
PLOT, xa, ya-y0, xr = [-1, 1], /xs, /ys, xtit = 'x', ytit = 'y', tit = 'y - S', thick = 2, $
      charsi = 1.4, yr = [-0.05, 0.05]
OPLOT, x, FLTARR(nx), psym = 1, symsi = 1.4
WRITE_PNG, 'example_CSN_1'+name+'.png', tvrd(/true)


WINDOW, 1, xsi = 500, ysi = 400
PLOT, x0, dy0, xr = [-1., 1], /xs, /ys, xtit = 'x', ytit = 'y', tit = 'dS/dx', thick = 2, $
      yrange = [-1.6, 2.2], charsi = 1.4
OPLOT, x[0:nx-1], b, psym = 1, symsi = 1.4
WRITE_PNG, 'example_CSN_2'+name+'.png', tvrd(/true)

WINDOW, 2, xsi = 500, ysi = 400
PLOT, x0, d2y0, xr = [-1., 1], /xs, /ys, xtit = 'x', ytit = 'y', tit = 'd!U2!NS/dx!U2!N', thick = 2, $
      yr = [-15, 15], charsi = 1.4
OPLOT, x[0:nx-1], 2*c, psym = 1, symsi = 1.4
WRITE_PNG, 'example_CSN_3'+name+'.png', tvrd(/true)

WINDOW, 3, xsi = 500, ysi = 400
PLOT, xa, ya, xr = [-2., 2], /xs, /ys, xtit = 'x', ytit = 'y', tit = 'S', thick = 2, $
      charsi = 1.4, yr = [-1, 1]


nxl = 201
xl = DINDGEN(nxl)/100-1.
cols = ['red', 'orange', 'yellow', 'green', 'blue', 'turquoise', 'magenta', 'pink']
; cls = GETCOLOR()
; cols = cls[FIX(randomu(seed, 50)*50)]
yl = DBLARR(nxl)
FOR j = 0, nx-2 DO BEGIN
  FOR i = 0, nxl-1 DO BEGIN  
    yl[i] = a[j] + b[j]*(xl[i]-x[j]) + c[j]*(xl[i]-x[j])^2 + d[j]*(xl[i]-x[j])^3
  ENDFOR
;   OPLOT, xl, yl, col = cols[j], lines = 2
  OPLOT, xl, yl, col = GETCOLOR(cols[j]), lines = 2
ENDFOR
OPLOT, x, y, psym = 1, symsi = 1.4
OPLOT, xa, ya
WRITE_PNG, 'example_CSN_4'+name+'.png', tvrd(/true)


; NO FRAME case
WINDOW, 3, xsi = 400, ysi = 200
PLOT, xa, ya, xr = [-1.2, 1.2], /xs, /ys,  $
      yr = [-1, 1], color = !p.background, $
      position = [0, 0, 1, 1]

nxl = 201
xl = DINDGEN(nxl)/100-1.
cols = ['red', 'orange', 'yellow', 'green', 'blue', 'turquoise', 'magenta', 'pink']
; cls = GETCOLOR()
; cols = cls[FIX(randomu(seed, 50)*50)]
yl = DBLARR(nxl)
FOR j = 0, nx-2 DO BEGIN
  FOR i = 0, nxl-1 DO BEGIN  
    yl[i] = a[j] + b[j]*(xl[i]-x[j]) + c[j]*(xl[i]-x[j])^2 + d[j]*(xl[i]-x[j])^3
  ENDFOR
;   OPLOT, xl, yl, col = cols[j], lines = 2
  OPLOT, xl, yl, col = GETCOLOR(cols[j]), lines = 2, thick = 2
ENDFOR
OPLOT, x, y, psym = 1, symsi = 1.7
OPLOT, xa, ya, thick = 2
WRITE_PNG, 'example_noborder.png', tvrd(/true)





stop
END
