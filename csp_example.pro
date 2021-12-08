pro csp_example

     ; The exact function for comparison.
     nxa = 241
     period = 2.
     xa = DINDGEN(nxa)/100-1.2
     ya = SIN(2.*!pi*xa/period)  
     

     ; The 'measurement' samples the function at 9 points (the array must be
     ; monotonically increasing).
     x = DOUBLE([-1., -0.8, -0.6, -0.45,  0.0, 0.3, 0.5, 0.6, 1.])
     nx = N_ELEMENTS(x)
     
     ; Force the exact periodicity, y[0] = y[nx-1]
     y = [SIN(2*!pi*x[0:nx-2]/period), SIN(2*!pi*x[0]/period)]
       
     ; Compute the coefficients of the approximative function S and evaluate S 
     ; at the same grid used to compute the exact function.
     nx0 = nxa
     cp = CSP_INTERPOLATION(x, y)
     x0 = FINDGEN(120)/4.-3
     y0 = CSP_EVALUATE(x0, x, cp)     
     
!p.background = GETCOLOR('white')
!p.color = GETCOLOR('black')

name = ''
; 
; ; Fig.1 The exact "unknown" function is shown in black. The locations of
; ;       the "measurements" are labeled with pluses. The interpolation
; ;       function is shown in red. 
WINDOW, 0, xsi = 500, ysi = 400

; The exact function for comparison.
nxm = 1001
period = 2.
xm = DINDGEN(nxm)/100-3.
ym = SIN(2.*!pi*xm/period)  

PLOT, xm, ym, xr = [-3, 7.], /xs, /ys, xtit = 'x', ytit = 'y', tit = 'S', thick = 2, $
      charsi = 1.4, yr = [-1.1, 1.1]
OPLOT, x, y, psym = 1, symsi = 2.
OPLOT, x0, y0, col = getcolor('red'), thick = 2
OPLOT, [-1, -1], [-1, 1], col = getcolor('dark gray'), lines = 2
OPLOT, [1, 1], [-1, 1], col = getcolor('dark gray'), lines = 2
WRITE_PNG, 'example_CSP_0'+name+'.png', tvrd(/true)




; Let's do the same but using the natural splines:
     nx0 = nxa
     cp = CSN_INTERPOLATION(x, y)
     x0 = FINDGEN(120)/4.-3
     y0 = CSN_EVALUATE(x0, x, y, cp) 

stop
END
