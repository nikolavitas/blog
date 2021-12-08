FUNCTION cs_function, x

; f = (SIN(0.3*x)+SIN(0.7*x)+SIN(1.1*x))/3
; f = x * SIN(2*!pi*x + 1)
f = 0.5 * x * COS(1.5*!pi*x + 0.5)
; f = x^3 + 2*x^2 + x + 3.

; f = 3*x + 4

; f = x^3.D0 * SIN(4.D0*!pi*x - 2.D0) + 0.8D0 * x* SIN(3.D0*!pi*x + 2.D0)

; f1 = SIN(3.*!pi*x)
; f2 = x * SIN(7.*!pi*x)
; f3 = 0.5*SIN(5.*!pi*x)
; f = f1 + f2 + f3

; Function that was example for the csn_interpolation application
; f = 0.5 * x * COS(1.5*!pi*x + 0.5)

; Function that was example for the csp_interpolation application
f = SIN(x)


RETURN, f
END
