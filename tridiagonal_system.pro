FUNCTION threediagonal_system, a, b, c, d
print, b
n = N_ELEMENTS(d)
x = DBLARR(n)

dprim = d
bprim = b

; My way
FOR i = 1, n-1 DO BEGIN
  dprim[i] = d[i] - (a[i]*d[i-1])/b[i-1]
  bprim[i] = b[i] - (a[i]*c[i-1])/b[i-1]
ENDFOR

x[n-1] = dprim[n-1]/bprim[n-1]
FOR i = n-2, 0, -1 DO x[i] = (dprim[i] - c[i]*x[i+1])/bprim[i]


; Wiki way
; c[0] = c[0]/b[0]
; d[0] = d[0]/b[0]
; 
; FOR i = 1, n-1 DO BEGIN
;   c[i] = c[i]/(b[i] - c[i-1]*a[i])
;   d[i] = (d[i] - d[i-1]*a[i])/(b[i] - c[i-1]*a[i])
; ENDFOR
; 
; x[n-1] = d[n-1]
; FOR i = n-2, 0, -1 DO x[i] = d[i] - c[i]*x[i+1]



RETURN, x
END