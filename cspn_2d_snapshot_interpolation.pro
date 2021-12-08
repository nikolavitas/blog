FUNCTION cspn_2d_snapshot_interpolation, x, y, snapshot

ns = SIZE(snapshot)
nx = ns[1]
ny = ns[2]

; Cubic spline interpolation
a_y = FLTARR(nx, ny-1) 
b_y = FLTARR(nx, ny-1) 
c_y = FLTARR(nx, ny-1) 
d_y = FLTARR(nx, ny-1) 

a_x = FLTARR(nx, ny) 
b_x = FLTARR(nx, ny) 
c_x = FLTARR(nx, ny) 
d_x = FLTARR(nx, ny) 

FOR ix = 0, nx-1 DO BEGIN
  column = REFORM(snapshot[ix, *])
  cv = CSN_INTERPOLATION(y, column)
  a_y[ix, *] = cv.a
  b_y[ix, *] = cv.b
  c_y[ix, *] = cv.c
  d_y[ix, *] = cv.d
ENDFOR

FOR iy = 0, ny-1 DO BEGIN
  row = REFORM(snapshot[*, iy])
  ch = CSP_INTERPOLATION([x, [x[0]]], [row, [row[0]]])
  a_x[*, iy] = ch.a
  b_x[*, iy] = ch.b
  c_x[*, iy] = ch.c
  d_x[*, iy] = ch.d
ENDFOR

result = {a_x:a_x, b_x:b_x, c_x:c_x, d_x:d_x, a_y:a_y, b_y:b_y, c_y:c_y, d_y:d_y}

RETURN, result

END
