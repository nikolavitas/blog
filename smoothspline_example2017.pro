PRO smoothspline_example2017
; This file is started in Feb 2017.
; I test splinecoeff, version that is public.
; I fixed a bug there that was noticed by someone in the user group.

nx = 100
x = findgen(100)                     
y = sin(x/2./!pi)*cos(x/2./!pi)^2    
p =  GAUSSIAN( findgen(nx), [1, 50, 10]) 
p = p/max(p)
z = splinecoeff(x, y, sigm = p, lambda = 1.d3)



; Plot the original data and the spline function for given lambda
; nx1 = 500
; x1 = FINDGEN(nx1)/100.+5.
; y1 = x1 * 0.

; nx1 = nx
; x1 = x
; y1 = y
; 
; FOR i = 0, nx1-1 DO begin
;   id = VALUE_LOCATE(x, x1[i])
;   y1[i] = z.d[id] + $
;           z.c[id] * (x1[i]-x[id]) + $
;           z.b[id] * (x1[i]-x[id])^2 + $
;           z.a[id] * (x1[i]-x[id])^3

; ENDFOR          

; ??
; y2 = csn_evaluate(x1, x, y, z)

y1 = spline_evaluate(x, x, z)
PLOT, x, y, col = 982982
OPLOT, x, y1
; OPLOT, x1, y2, psym = 2
stop

END
