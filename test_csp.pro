pro test_csp

     ; The exact function for comparison.
     nxa = 201
     xa = DINDGEN(nxa)/100-1.
     ya = SIN(!pi*xa)  

     ; The 'measurement' samples the function at 9 points (the array must be
     ; monotonically increasing).
     x = DOUBLE([-1., -0.8, -0.6, -0.45,  0.0, 0.3, 0.5, 0.6, 1.])
     nx = N_ELEMENTS(x)
     
     ; Force the exact periodicity, y[0] = y[nx-1]
     y = [SIN(!pi*x[0:nx-2]), SIN(!pi*x[0])]
       
     ; Compute the coefficients of the approximative function S and evaluate S 
     ; at the same grid used to compute the exact function.
     x0 = xa
     nx0 = nxa
     cp = CSP_INTERPOLATION(x, y)
     
     y0 = DBLARR(nx0)
     FOR i = 0, nx0-1 DO BEGIN
       j = VALUE_LOCATE(x, x0[i]) 
       IF (j GE 0) AND (j LT nx-1) THEN BEGIN
         y0[i] = cp.a[j] + cp.b[j]*(x0[i]-x[j]) + cp.c[j]*(x0[i]-x[j])^2 + cp.d[j]*(x0[i]-x[j])^3
       ENDIF 
     ENDFOR

     plot, xa, ya
     oplot, x, y, psym = 1
     oplot, x0, y0
     
; 
; x = FINDGEN(7)
; y = SIN(x)
; x = [x, 2*!pi]
; y = [y, 0]
; 
; co = csp_interpolation(x, y)
; ; cn = csn_interpolation(x, y)
; 
; 
; xx = FINDGEN(61)/10.
; yyp = FLTARR(61)
; yyn = FLTARR(61)
; for i = 0, 60 do begin
;   j = VALUE_LOCATE(x, xx[i])
;   print, i, xx[i], j, x[j]
;   yyp[i] = co.a[j]+ co.b[j]*(xx[i]-x[j])  + co.c[j]*(xx[i]-x[j])^2  + co.d[j]*(xx[i]-x[j])^3
; ;   yyn[i] = cn.a[j]+ cn.b[j]*(xx[i]-x[j])  + cn.c[j]*(xx[i]-x[j])^2  + cn.d[j]*(xx[i]-x[j])^3
; endfor
; ; 
; r = IMSL_csinterp(x, y, /periodic)
; 
; a = r.coef[0:27:4]
; b = r.coef[1:27:4]
; c = r.coef[2:27:4]
; d = r.coef[3:27:4]
; 
; ; yyi = FLTARR(61)
; ; yy2 = FLTARR(61)
; ; yy3 = FLTARR(61)
; ; yy1 = FLTARR(61)
; ; 
; ; for i = 0, 60 do begin
; ;   j = VALUE_LOCATE(x, xx[i])
; ;   print, i, xx[i], j, x[j]
; ;   yyi[i] = a[j]+ b[j]*(xx[i]-x[j])  + c[j]*(xx[i]-x[j])^2  + d[j]*(xx[i]-x[j])^3
; ;   yy1[i] = b[j]  + 2*c[j]*(xx[i]-x[j])  + 3*d[j]*(xx[i]-x[j])^2
; ;   yy2[i] = 2*c[j]  + 6*d[j]*(xx[i]-x[j]) 
; ;   yy3[i] = 6*d[j]
; ; endfor
; ; 
; yyy = IMSL_SPVALUE(xx, r)
; 
; plot, x, y, psym = 1
; oplot, xx, yyp
; ; oplot, xx, yyn, lines = 0, col = getcolor('red')
; oplot, xx, yyy, lines = 1, col = getcolor('green')

stop
END
