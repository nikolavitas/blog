FUNCTION trace_rays_3d, ix0 = ix0, theta = theta, threshold = threshold
;theta = 85.
;ix0 = 0
;threshold = 1.; 49.99999

; If a section point with a vertical wall is closer to a
; section with horizontal wall than the threshold (relative to 
; the spacing), than this section point is disregarded.

; The routine apparently misbehaves for small threshold values, thus when 
; the sections with the vertical walls are included. 

;-------------------------------------------------------------------------------
; Set the axes of the snapshot
;-------------------------------------------------------------------------------
ny = 100
dy = 14.
;dy = 20.
y = FINDGEN(ny) * dy


nx = 288
dx = 6000./nx
;dx = 20.
x = FINDGEN(nx) * dx

nz = nx

;-------------------------------------------------------------------------------
; Load the snapshot
;-------------------------------------------------------------------------------
dir = '/Users/nikola/Data/000000/dat_mu1.00/'
t3d = FLTARR(nx, ny, nz)
OPENR, 1, dir + 't3d.000000.dat'
READU, 1, t3d
CLOSE, 1

;-------------------------------------------------------------------------------
; Set the region at the upper boundary from which the ray are traced.
;-------------------------------------------------------------------------------
xr = [0, 287]
zr = [0, 287]

;-------------------------------------------------------------------------------
; Interpolation along the grid using the bilinear method
;-------------------------------------------------------------------------------
; Here I don't use trace_ray function. The interpolation points are selected
; equidistantly along the ray, no matter if they pass through a grid cell
; or not. The input for tracing is y0 (the horizontal snapshot at which
; I start the rays (in both directions if y0 NE 0 or y0 NE ny-1)), theta
; and ds along the ray.

method = 'bilinear'

IF method EQ 'bilinear' THEN BEGIN

  ds = 1. ; km
  dsx = ds * SIN(theta*!pi/180.) ; change along x-axis while s along the path changes for ds
  dsy = ds * COS(theta*!pi/180.) ; change along y-axis while s along the path changes for ds

  ; Number of points along the ray
  ymax = MAX(y) 
  nd = ymax/dsy

  ; s is length along the ray
  s = DINDGEN(nd)*ds

  t3d_mu = REFORM(DBLARR(xr[1]-xr[0]+1, nd, zr[1]-zr[0]+1), xr[1]-xr[0]+1, nd, zr[1]-zr[0]+1)

  ; sx and sy are projections of s onto x and y axes.
  sx = s * SIN(theta*!pi/180.)
  sy = s * COS(theta*!pi/180.)

  ; xc and yc are the coordinates of grid cells that contain the 
  ; points along the ray
  xc = FLOOR(sx / dx)
  yc = FLOOR(sy / dy)

  ; every time the x coordinate goes over nx we lower it for nx, 
  ; so to keep it in the interval between 0 and nx-1).
  FOR id = 0, nd-1 DO $
    WHILE xc[id] GE nx DO xc[id] = xc[id] - nx
    
  FOR iz0 = zr[0], zr[1] DO BEGIN 

    ; Extract a 2D slice that will be interpolated
    t2d = REFORM(t3d[*, *, iz0])

    ng = 1 ; number of ghost cells
 
    ; we extend the 2D slice and the x-axis for one ghost cell/row 
    ; (enough to accomodate the bilinear interpolation).  
    t2d_extended = DBLARR(nx+ng, ny)
    t2d_extended[0:nx-1, *] = t2d
    t2d_extended[nx-1+ng, *] = t2d[0, *]
    x_extended = [x, x[nx-1]+dx]

    FOR ix0 = xr[0], xr[1] DO BEGIN
      FOR id = 0, nd-1 DO BEGIN

        ; xc is valid only for the rays starting from x = 0. For 
        ; other rays we add ix0 and call it xc0
        xc0 = xc[id] + ix0

        nperiod = 0
        WHILE xc0 GE nx DO BEGIN
          xc0 = xc0 - nx
          nperiod += 1
        ENDWHILE

        zz = t2d_extended[xc0:xc0+1, yc[id]:yc[id]+1]
       ; print, zz[0, 0], zz[1, 1]    


; This is fine except in the corner
        
        x0 = LAGRANGE(sx[id]+(ix0-nperiod*nx)*dx, REFORM(zz[*, 0]), x_extended[xc0:xc0+1])
        x1 = LAGRANGE(sx[id]+(ix0-nperiod*nx)*dx, REFORM(zz[*, 1]), x_extended[xc0:xc0+1])

        t3d_mu[ix0-xr[0], id, iz0-zr[0]] = LAGRANGE(sy[id], [x0, x1], y[yc[id]:yc[id]+1])
    
      ENDFOR
    ENDFOR
  ENDFOR
ENDIF
;-------------------------------------------------------------------------------

stop

;-------------------------------------------------------------------------------
; Evaluation of the splines at the sections of the ray and the grid
;-------------------------------------------------------------------------------
IF method EQ 'splines' THEN BEGIN

;-------------------------------------------------------------------------------
; Trace a tray starting from an "upper left" grid point and
; propagating "down" and "right". Ray tracing does not depend on the content,
; but only on the grid size and spacing.
;-------------------------------------------------------------------------------
path = TRACE_RAY(dx, dy, nx, ny, theta, /remove, threshold = threshold)

nd = N_ELEMENTS(path.xs)   ; Number of points along the ray
wall = path.wall           ; Indicator if the wall is (v)ertical or (h)orizontal
xs = path.xs               ; Length of the ray path along x
ys = path.ys               ; Length of the ray path along y
xc = path.xc               ; Coordinate of the wall that is crossed in x
yc = path.yc               ; Coordinate of the wall that is crossed in y

;-------------------------------------------------------------------------------
; Define the new array that will contain the skewed snapshot
;-------------------------------------------------------------------------------
t3d_mu = REFORM(DBLARR(xr[1]-xr[0]+1, nd, zr[1]-zr[0]+1), xr[1]-xr[0]+1, nd, zr[1]-zr[0]+1)


FOR iz0 = zr[0], zr[1] DO BEGIN 

  ; Extract a 2D slice that will be interpolated
  t2d = REFORM(t3d[*, *, iz0])

  ;-------------------------------------------------------------------------------
  ; Cubic spline computation
  ; Splines depend on the content and on the grid, but not on the ray tracing.
  ; Here we compute only the coefficient! We evaluate the splines in the next
  ; step. 
  ;-------------------------------------------------------------------------------      
  coeff = CSPN_2D_SNAPSHOT_INTERPOLATION(x, y, t2d)

  a_x = coeff.a_x
  b_x = coeff.b_x
  c_x = coeff.c_x
  d_x = coeff.d_x

  a_y = coeff.a_y
  b_y = coeff.b_y
  c_y = coeff.c_y
  d_y = coeff.d_y

  FOR ix0 = xr[0], xr[1] DO BEGIN
    FOR id = 0, nd-1 DO BEGIN

      ; How many times the ray passed the entire domain?
      nperiod = FIX((xs[id] + ix0*dx)/(nx*dx))
   
      ; The coordinate of the current cell in which the ray crosses the grid.
      xc0 = xc[id] - nperiod*nx + ix0

      ; If horizontal wall is crossed, do the spline evaluation using the
      ; coefficient _x computed for the corresponding horizontal grid line
      IF wall[id] EQ 'h' THEN BEGIN 
        xs0 = xs[id]+ix0*dx-nperiod*nx*dx
        t3d_mu[ix0-xr[0], id, iz0-zr[0]] = $
              a_x[xc0, yc[id]] + $
              b_x[xc0, yc[id]]*(xs0-x[xc0]) + $
              c_x[xc0, yc[id]]*(xs0-x[xc0])^2 + $
              d_x[xc0, yc[id]]*(xs0-x[xc0])^3 
      ENDIF

      ; If vertical wall is crossed, do the spline evaluation using the
      ; coefficient _y computed for the corresponding vertical grid line
      IF wall[id] EQ 'v' THEN BEGIN
        ys0 = (ny-1)*dy - ys[id] ; distance starting from the bottom
        t3d_mu[ix0-xr[0], id, iz0-zr[0]] = $
              a_y[xc0, yc[id]] + $
              b_y[xc0, yc[id]]*(ys0-y[yc[id]]) + $
              c_y[xc0, yc[id]]*(ys0-y[yc[id]])^2 + $
              d_y[xc0, yc[id]]*(ys0-y[yc[id]])^3 
      ENDIF

    ENDFOR
  ENDFOR
ENDFOR
ENDIF
;-------------------------------------------------------------------------------


; The output contains the interpolated quantity, the length along the ray (s),
; and the x axis that remains unchanged. 

a = {quantity:t3d_mu, s:SQRT(xs^2 + ys^2), x:x}

RETURN, a

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
