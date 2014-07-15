pro jeje_symbols,nsym,scale,color=color

; symbols are somewhat normalised to look the same size... (side ~ 1) 

; 1/2 : open/filled circle
; 3/4 : open/filled square
; 5/6 : open/filled downwards triangle
; 7/8 : open/filled upwards triangle

fill = 0

case nsym of

    1: begin   ; open circle
        a = findgen(30)
        a = a * (2.*!pi / 28.)
        xarr = cos(a) / 2.0 
        yarr = sin(a) / 2.0
    end
 
    2: begin   ; filled circle
        a = findgen(30)
        a = a * (2.*!pi / 28.)
        xarr = cos(a) / 2.0
        yarr = sin(a) / 2.0
        fill = 1
    end

    3: begin  ; open square box
        xarr = ( yarr = fltarr(5)) 
        xarr(0) = -1 & yarr(0) = -1 
        xarr(1) = 1  & yarr(1) = -1 
        xarr(2) = 1  & yarr(2) = 1
        xarr(3) = -1 & yarr(3) = 1
        xarr(4) = -1 & yarr(4) = -1 
        xarr = xarr / 2.0
        yarr = yarr / 2.
    end
    
    4: begin  ; filled square box
        xarr = ( yarr = fltarr(5)) 
        xarr(0) = -1 & yarr(0) = -1 
        xarr(1) = 1  & yarr(1) = -1 
        xarr(2) = 1  & yarr(2) = 1
        xarr(3) = -1 & yarr(3) = 1
        xarr(4) = -1 & yarr(4) = -1 
        xarr = xarr / 2.
        yarr = yarr / 2.
        fill = 1
    end

    5: begin  ; open downwards triangle
       xarr = ( yarr = fltarr(4))
       xarr(0) = 0.     & yarr(0) = 1.0
       xarr(1) = 0.866  & yarr(1) = - 0.5
       xarr(2) = -0.866 & yarr(2) = - 0.5
       xarr(3) = 0.0    & yarr(3) = 1.0
       xarr = xarr / 1.4
       yarr = (yarr-0.25) / 1.4
    end

    6: begin  ; filled downwards triangle
       xarr = ( yarr = fltarr(4))
       xarr(0) = 0.     & yarr(0) = 1.0
       xarr(1) = 0.866  & yarr(1) = - 0.5
       xarr(2) = -0.866 & yarr(2) = - 0.5
       xarr(3) = 0.0    & yarr(3) = 1.0
       xarr = xarr / 1.4
       yarr = (yarr-0.25) / 1.4
       fill = 1
    end

    7: begin  ; open upwards triangle
       xarr = ( yarr = fltarr(4))
       xarr(0) = 0.     & yarr(0) = -1.0
       xarr(1) = 0.866  & yarr(1) = 0.5
       xarr(2) = -0.866 & yarr(2) = 0.5
       xarr(3) = 0.0    & yarr(3) = -1.0
       xarr = xarr / 1.4
       yarr = (yarr-0.25) / 1.4
    end

    8: begin  ; filled upwards triangle
       xarr = ( yarr = fltarr(4))
       xarr(0) = 0.     & yarr(0) = -1.0
       xarr(1) = 0.866  & yarr(1) = 0.5
       xarr(2) = -0.866 & yarr(2) = 0.5
       xarr(3) = 0.0    & yarr(3) = -1.0
       xarr = xarr / 1.4
       yarr = (yarr-0.25) / 1.4
       fill = 1
    end
endcase

xarr = xarr * scale
yarr = yarr * scale

;set symbol buffer
usersym,xarr,yarr,fill=fill,color=color,thick=3

return
end
