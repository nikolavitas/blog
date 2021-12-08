FUNCTION trace_ray, dx, dy, nx, ny, theta, remove = remove, threshold = threshold
thetapi = theta*!pi/180.

; dsx = Path that ray travels along the x-direction corresponding to the grid-cell height dy
; dsy = Path that ray travels along the y-direction corresponding to the grid-cell width
dsy = dx/TAN(thetapi)
dsx = TAN(thetapi)*dy

x = DINDGEN(nx)*dx
y = DINDGEN(ny)*dy

; We determine the horizontal crossing points (just from the fact that we know how much the ray 
; has to travel between them).
sx0 = (DINDGEN(ny))*dsx

; Previous coordinate of the cell in which a horizontal section took place.
ix0 = FIX(sx0/dx)

; The corresponding y values are simply 
sy0 = (DINDGEN(ny))*dy

;FOR iy = 0, ny-1 DO $
;  OPLOT, [sx0[iy]], [(ny-1)*dy - sy0[iy]], col = GETCOLOR('turquoise'), psym = -1, symsi = 3

; Initiate the arrays
xs = [sx0[0]]
ys = [sy0[0]]
xc = [0]
yc = [ny-1]
wall = ['h']

FOR iy = 1, ny-1 DO BEGIN

  ; nvsect is the number of vertical wall crossed from the last crossing of 
  ; a horizontal wall.
  nvsect = ix0[iy] - ix0[iy-1]
  
  ; If the number of vertical sections is larger than zero, add them all to the 
  ; list of section points.
  IF nvsect GT 0 THEN BEGIN

    ; The indices of these walls are:
    ix1 = INDGEN(nvsect) + ix0[iy-1] + 1 

    ; Their x-coordinates are:
    sx1 = ix1*dx

    ; Their y-coordinates are:
    sy1 = sx1/TAN(theta*!pi/180.)

    ; Upgrade the arrays with the sections of the vertical walls
    xs = [xs, sx1]  
    ys = [ys, sy1]  
    xc = [xc, ix1]
    yc = [yc, ny-1-(INTARR(nvsect)+iy)]
    wall = [wall, STRARR(nvsect) + 'v']   

;    FOR iyy = 0, nvsect-1 DO BEGIN
;      OPLOT, [sx1[iyy]], [(ny-1)*dy - sy1[iyy]], col = GETCOLOR('lime green'), psym = -1, symsi = 3
;    ENDFOR

  ENDIF

  ; Add the sections of the horizontal walls to the arrays.
  xs = [xs, sx0[iy]]  
  ys = [ys, sy0[iy]] 
  xc = [xc, ix0[iy]]
  yc = [yc, ny-1-iy]
  wall = [wall, 'h']
  
ENDFOR


nd = N_ELEMENTS(xs)

; Remove the sections of the vertical walls that are closer than a given
; threshold to the nearest other point 
IF NOT(KEYWORD_SET(threshold)) THEN BEGIN
  threshold = 1. ; in % of the flight path if there is no vertical walls
ENDIF

IF KEYWORD_SET(remove) THEN BEGIN
  shortch = SQRT(dsx^2 + dy^2)
  id = 1
  WHILE id LT nd DO BEGIN
    distance = SQRT((xs[id]-xs[id-1])^2 + (ys[id]-ys[id-1])^2)
    IF distance/shortch*100 LT threshold THEN BEGIN
      IF wall[id] EQ 'h' THEN BEGIN
        ; If it's a section of a horizontal wall, then remove the
        ; previous point (that must be a section of the horizontal
        ; wall if threshold is LT 100.
        xs = [xs[0:id-2], xs[id:nd-1]]
        ys = [ys[0:id-2], ys[id:nd-1]]
        xc = [xc[0:id-2], xc[id:nd-1]]
        yc = [yc[0:id-2], yc[id:nd-1]]
        wall = [wall[0:id-2], wall[id:nd-1]]
      ENDIF ELSE BEGIN
        PRINT, 'remove [id]', id, wall[id]
        xs = [xs[0:id-1], xs[id+1:nd-1]]
        ys = [ys[0:id-1], ys[id+1:nd-1]]
        xc = [xc[0:id-1], xc[id+1:nd-1]]
        yc = [yc[0:id-1], yc[id+1:nd-1]]
        wall = [wall[0:id-1], wall[id+1:nd-1]]
      ENDELSE 
      id -= 1
      nd = N_ELEMENTS(xs)
    ENDIF
    id += 1
  ENDWHILE
ENDIF

s = DBLARR(nd-1)
FOR id = 0, nd-2 DO BEGIN
  s[id] = SQRT((xs[id+1]-xs[id])^2 + (ys[id+1]-ys[id])^2)
ENDFOR

path = {s:s, xs:xs, ys:ys, xc:xc, yc:yc, wall:wall}

RETURN, path
END
