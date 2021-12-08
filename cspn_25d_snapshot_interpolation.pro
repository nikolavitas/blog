FUNCTION cspn_25d_snapshot_interpolation, x, y, z, snapshot

ns = SIZE(snapshot)
nx = ns[1]
ny = ns[2]
nz = ns[3]

; Cubic spline interpolation
a_y = FLTARR(nx, nz, ny-1) 
b_y = FLTARR(nx, nz, ny-1) 
c_y = FLTARR(nx, nz, ny-1) 
d_y = FLTARR(nx, nz, ny-1) 

a_x = FLTARR(nx, nz, ny) 
b_x = FLTARR(nx, nz, ny) 
c_x = FLTARR(nx, nz, ny) 
d_x = FLTARR(nx, nz, ny) 


FOR ix = 0, nx-1 DO BEGIN
  FOR iz = 0, nz-1 DO BEGIN
    column = REFORM(snapshot[ix, *, iz])
    cy = CSN_INTERPOLATION(y, column)
    a_y[ix, iz, *] = cy.a
    b_y[ix, iz, *] = cy.b
    c_y[ix, iz, *] = cy.c
    d_y[ix, iz, *] = cy.d
  ENDFOR
ENDFOR

FOR iy = 0, ny-1 DO BEGIN
  FOR iz = 0, nz-1 DO BEGIN
    row = REFORM(snapshot[*, iy, iz])
    cx = CSP_INTERPOLATION([x, [x[0]]], [row, [row[0]]])
    a_x[*, iz, iy] = cx.a
    b_x[*, iz, iy] = cx.b
    c_x[*, iz, iy] = cx.c
    d_x[*, iz, iy] = cx.d
  ENDFOR
ENDFOR

result = {a_x:a_x, b_x:b_x, c_x:c_x, d_x:d_x, $
              a_y:a_y, b_y:b_y, c_y:c_y, d_y:d_y}

RETURN, result

END
