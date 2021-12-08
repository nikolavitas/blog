FUNCTION ray_tracing_example, ix0 = ix0, theta = theta, threshold = threshold
;theta = 85.
;ix0 = 0
;threshold = 1.; 49.99999

ny = 100
dy = 14.
;dy = 20.
y = FINDGEN(ny) * dy


nx = 288
dx = 6000./nx
;dx = 20.
x = FINDGEN(nx) * dx

nz = nx
dz = dx
z = x

dir = '/Users/nikola/Data/000000/dat_mu1.00/'
t3d = FLTARR(nx, ny, nz)
OPENR, 1, dir + 't3d.000000.dat'
READU, 1, t3d
CLOSE, 1

t2d = REFORM(t3d[*, *, 0])

xr = [0, 10]
zr = [0, 10]


;-------------------------------------------------------------------------------
; Trace a tray starting from an "upper left" grid point and
; propagating "down" and "right"
;-------------------------------------------------------------------------------
path = TRACE_RAY(dx, dy, nx, ny, theta, /remove, threshold = threshold)

nd = N_ELEMENTS(path.xs)
t1d_mu = DBLARR(nd)
wall = path.wall
xs = path.xs
ys = path.ys
xc = path.xc
yc = path.yc
;-------------------------------------------------------------------------------

;-------------------------------------------------------------------------------
; Cubic spline interpolation
;-------------------------------------------------------------------------------      
;coeff = CSPN_2D_SNAPSHOT_INTERPOLATION(x, y, t2d)
coeff = CSPN_25D_SNAPSHOT_INTERPOLATION(x, y, z, t3d)

a_x = coeff.a_x
b_x = coeff.b_x
c_x = coeff.c_x
d_x = coeff.d_x
a_y = coeff.a_y
b_y = coeff.b_y
c_y = coeff.c_y
d_y = coeff.d_y
;-------------------------------------------------------------------------------

;-------------------------------------------------------------------------------
; Evaluation of the splines at the sections of the ray and the grid
;-------------------------------------------------------------------------------
t1d_mu = DBLARR(xr[1]-xr[0]+1, zr[1]-zr[0]+1, nd)
FOR iz0 = zr[0], zr[1] DO BEGIN
FOR ix0 = xr[0], xr[1] DO BEGIN
  FOR id = 0, nd-1 DO BEGIN

    ; How many times the ray passed the entire domain?
    nperiod = FIX((xs[id] + ix0*dx)/(nx*dx))
   
    ; The coordinate of the current cell in which the ray crosses the grid.
    xc0 = xc[id] - nperiod*nx + ix0

    IF wall[id] EQ 'h' THEN BEGIN 
      xs0 = xs[id]+ix0*dx-nperiod*nx*dx
      t1d_mu[ix0, iz0, id] = a_x[xc0, iz0, yc[id]] + $
                             b_x[xc0, iz0, yc[id]]*(xs0-x[xc0]) + $
                             c_x[xc0, iz0, yc[id]]*(xs0-x[xc0])^2 + $
                             d_x[xc0, iz0, yc[id]]*(xs0-x[xc0])^3 
    ENDIF

    IF wall[id] EQ 'v' THEN BEGIN
      ys0 = (ny-1)*dy - ys[id] ; distance starting from the bottom
      t1d_mu[ix0, iz0, id] = a_y[xc0, iz0, yc[id]] + $
                            b_y[xc0, iz0, yc[id]]*(ys0-y[yc[id]]) + $
                            c_y[xc0, iz0, yc[id]]*(ys0-y[yc[id]])^2 + $
                            d_y[xc0, iz0, yc[id]]*(ys0-y[yc[id]])^3 
     ENDIF

  ENDFOR
ENDFOR
ENDFOR
;-------------------------------------------------------------------------------


;stop

a = {y:t1d_mu, x:xs}

RETURN, a

;stop
END
;nx1 = 10
;ny1 = 10
;!p.color = GETCOLOR('black')
;!p.background = GETCOLOR('white')
;WINDOW, 0, xsi = 600, ysi = 600
;PLOT, [0], [0], xr = [0, nx1*dx], yr = [ (ny-ny1-1)*dy, (ny-1)*dy], ys = 1, xs = 1
;FOR i = 0, nx1 DO OPLOT, [i*dx, i*dx], [ (ny-ny1-1)*dy, (ny-1)*dy]
;FOR j = (ny-ny1-1), ny-1 DO OPLOT, [0, nx*dx], [j*dy, j*dy]
;FOR iy = 0, ny-1 DO $
;  OPLOT, [path.xs[iy]], [(ny-1)*dy - path.ys[iy]], col = GETCOLOR('turquoise'), psym = -1, symsi = 3

; Now we do the interpolation along the grid lines at every
; intersection point.
