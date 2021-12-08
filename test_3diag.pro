PRO test_3diag

; x0 = FINDGEN(5)
; aa = DOUBLE([6., 2.,  3., 4., 1.])
; bb = DOUBLE([3., 4., 11., 7., 2.])
; cc = DOUBLE([1., 1.,  1., 3., 3.])
; 
; dd = DBLARR(5)
; dd[0] = bb[0]*x0[0] + cc[0]*x0[1] + aa[0]*x0[4]
; FOR i = 1, 3 DO dd[i] = aa[i]*x0[i-1] + bb[i]*x0[i] + cc[i]*x0[i+1]
; dd[4] = aa[4]*x0[3] + bb[4]*x0[4] + cc[4]*x0[0]
; 
; x = threediagonal_periodic_system(aa, bb, cc, dd)
; 
; print, x[0]*bb[0] + x[1]*cc[0] + x[4]*aa[0], dd[0]
; FOR i = 1, 3 DO PRINT, x[i-1]*aa[i] + x[i]*bb[i] + x[i+1]*cc[i], dd[i]
; print, x[3]*aa[4] + x[4]*bb[4] + x[0]*cc[4], dd[4]



x0 = FINDGEN(5)
aa = [0., 2.,  3., 4., 1.]
bb = [3., 4., 11., 7., 2.]
cc = [1., 1.,  1., 3., 0.]

dd = FLTARR(5)
dd[0] = bb[0]*x0[0] + cc[0]*x0[1] + aa[0]*x0[4]
FOR i = 1, 3 DO dd[i] = aa[i]*x0[i-1] + bb[i]*x0[i] + cc[i]*x0[i+1]
dd[4] = aa[4]*x0[3] + bb[4]*x0[4] + cc[4]*x0[0]

x = threediagonal_system(aa, bb, cc, dd)

print, x

print, x[0]*bb[0] + x[1]*cc[0] + x[4]*aa[0], dd[0]
FOR i = 1, 3 DO PRINT, x[i-1]*aa[i] + x[i]*bb[i] + x[i+1]*cc[i], dd[i]
print, x[3]*aa[4] + x[4]*bb[4] + x[0]*cc[4], dd[4]


END


; ; We form a system of five equations
; aa = [0., 2.,  3., 4., 1.]
; bb = [3., 4., 11., 7., 2.]
; cc = [1., 1.,  1., 3., 0.]
; 
; ; Set the exact solution
; x0 = FINDGEN(5)
; 
; ; Compute the right-hand side to complete the system
; dd = FLTARR(5)
; dd[0] = bb[0]*x0[0] + cc[0]*x0[1] + aa[0]*x0[4]
; FOR i = 1, 3 DO dd[i] = aa[i]*x0[i-1] + bb[i]*x0[i] + cc[i]*x0[i+1]
; dd[4] = aa[4]*x0[3] + bb[4]*x0[4] + cc[4]*x0[0]
; 
; ; Solve the system
; x = threediagonal_system(aa, bb, cc, dd)
; 
; ; Compare the solution of the system and the exact solution
; PRINT, x
; PRINT, x0



; ; We form a system of five equations
; IDL> aa = [6., 2.,  3., 4., 1.]
; IDL> bb = [3., 4., 11., 7., 2.]
; IDL> cc = [1., 1.,  1., 3., 3.]
; 
; ; Set the exact solution
; IDL> x0 = FINDGEN(5)
; 
; ; Compute the right-hand side to complete the system
; IDL> dd = FLTARR(5)
; IDL> dd[0] = bb[0]*x0[0] + cc[0]*x0[1] + aa[0]*x0[4]
; IDL> FOR i = 1, 3 DO dd[i] = aa[i]*x0[i-1] + bb[i]*x0[i] + cc[i]*x0[i+1]
; IDL> dd[4] = aa[4]*x0[3] + bb[4]*x0[4] + cc[4]*x0[0]
; 
; ; Solve the system
; IDL> x = threediagonal_periodic_system(aa, bb, cc, dd)
; 
; ; Compare the solution of the system and the exact solution
; IDL> PRINT, x
;       0.00000      1.00000      2.00000      3.00000      4.00000
; IDL> PRINT, x0
;       0.00000      1.00000      2.00000      3.00000      4.00000

