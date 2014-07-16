pro plot_all_list,dir,method,feps
  
  winnb = 0

  create_tree_list,dir,tree_out
  
  for itree = 0,n_elements(tree_out) - 1 do begin
     plot_full_merger_tree2,dir,tree_out(itree),method,feps,winnb
  end
  
end

;---------------------------------------------------------------------------------
pro create_tree_list,dir,tree_out
  
  print,"create_tree_list"
  
  read_where_haloes_are,dir,tree_numb,halo_numb,st_numb,branch_numb,branch_id,st_start,st_end

  n_tree_out = 0
  for itree = 0,n_elements(tree_numb) - 1 do begin
     found = 0
     if(itree gt 0) then begin
        for jtree = 0,itree-1 do begin
           if(tree_numb(jtree) eq tree_numb(itree)) then found = 1
        end
     endif
     if(found eq 0) then n_tree_out = n_tree_out + 1
  end
  
  tree_out = intarr(n_tree_out)
  n_tree_out = 0
  for itree = 0,n_elements(tree_numb) - 1 do begin
     found = 0
     if(itree gt 0) then begin
        for jtree = 0,itree-1 do begin
           if(tree_numb(jtree) eq tree_numb(itree)) then found = 1
        end
     endif
     if(found eq 0) then begin
        tree_out(n_tree_out) = tree_numb(itree)
        n_tree_out           = n_tree_out + 1
     end
  end
  
  print,'number of full trees to show:',n_tree_out

end

;---------------------------------------------------------------------------------
pro plot_full_merger_tree2,dir,tree_number,method,feps,winnb

  file_tree = string(dir,"full_tree_",tree_number,"_",method,".dat",format='(a,a,i8.8,a1,a3,a)')
  if(file_test(file_tree) ne 1) then begin
     print,"file ",file_tree," doesn't exist"
     return
  endif 
  read_tree, file_tree,br_id,st_id,son,host,mass,f_lev,f_frag
  
  ibmax   = max(br_id)
  ibmarg  = ibmax/20 + 1
  ibmax   = ibmax + ibmarg
  ibmin   = 1 - ibmarg
  
  stmax   = max(st_id)
  st_id(0)= stmax
  stmin   = min(st_id)
  st_id(0)= 0
  st_up   = stmax + (stmax-stmin)/20. 
  dst_leg = (stmax-stmin)/20.

  xplotsize  = 1.5

  if(ibmax ge 250) then xplotsize = 2.
  if(ibmax ge 500) then xplotsize = 2.5
  if(ibmax ge 750) then xplotsize = 3.
  if(ibmax ge 1000) then xplotsize = 3.5
  if(ibmax ge 1500) then xplotsize = 4.
  if(ibmax ge 2500) then xplotsize = 4.5
  if(ibmax ge 5000) then xplotsize = 5.
  
  if(stmax le 50) then begin
     yplotsize = 1.
  endif else yplotsize = 1.5
  

  col = [6,7,12]
  loadct,12,ncolors=16
  file_eps = string(dir,"full_tree_",tree_number,"_",method,".eps",format='(a,a,i8.8,a1,a3,a)')
  ; xsize,ysize,!p.charsize,!x.charsize,!y.charsize,!p.charthick,feps,window,file
  init_plot,xplotsize,yplotsize,1.5,1.5,1.5,5.,feps,winnb,file_eps
  plot,[ibmin,ibmax],[stmin,stmax],xstyle=1,ystyle=1,/nodata,xmargin=[5,5],ymargin=[5,5]
  for i = 1l,n_elements(br_id)-1 do begin
     ison = son(i)
     if(ison gt 0) then begin
        icol = 0
        if(f_lev(i) gt 1)  then icol = 1
       ; if(f_frag(i) gt 0) then icol = 2
        oplot,[br_id(i),br_id(ison)],[st_id(i),st_id(ison)],col=col(icol)
     endif else begin
        ihost = host(i)
        if(ihost gt 0) then begin
           plots,[br_id(i),br_id(i)],[st_id(i),st_up]
           plots,[br_id(ihost),br_id(i)],[st_id(ihost),st_up]
        end
     endelse
     if(f_frag(i) eq 1) then oplot,[br_id(i)],[st_id(i)],psym=7
  end

  find_targets,dir,tree_number,n_targ,halo_targnumb,st_targnumb,br_targnumb
  if(n_targ gt 0) then begin
     for itarg = 0 ,n_targ - 1 do begin
        point_target,br_targnumb(itarg),st_targnumb(itarg),dst_leg,halo_targnumb(itarg),col(2)
     end
  end
  
  if(feps eq 1) then begin
     device,/close
     set_plot,"x"
  endif
  
  if(n_targ gt 0) then begin
     for itarg = 0 ,n_targ - 1 do begin
        plot_zoom_merger_tree2,dir,tree_number,halo_targnumb(itarg),st_targnumb(itarg),method,feps,winnb
     end
  end
  
end

;---------------------------------------------------------------------------------
pro plot_zoom_merger_tree2,dir,tree_number,halo_number,st_number,method,feps,winnb
 
  read_where_haloes_are,dir,tree_numb,halo_numb,st_numb,branch_numb,branch_id,st_start,st_end
  

  col = [6,7,12]

  for itree = 0,n_elements(tree_numb) - 1 do begin
     if(tree_numb(itree) eq tree_number and halo_numb(itree) eq halo_number and st_numb(itree) eq st_number) then begin
        stz_start = st_start(itree)
        stz_end   = st_end(itree)
        br_number = branch_id(itree)

        file_zoom = string(dir,"zoom_tree_",tree_number,"_",halo_number,"_",st_number,"_",stz_start,"_",stz_end,"_",method,".dat",format='(a,a,2(i8.8,a1),3(i3.3,a1),a3,a)')
        
        read_tree, file_zoom,br_id,st_id,son,host,mass,f_lev,f_frag
        

        ibmax   = max(br_id)
        ibmarg  = ibmax/20 + 1
        ibmax   = ibmax + ibmarg
        ibmin   = 1 - ibmarg
        stmax   = max(st_id)
        st_id(0)= stmax
        stmin   = min(st_id)
        st_id(0)= 0
        st_up   = stmax + (stmax-stmin)/20. 
        dst_leg = (stmax-stmin)/20.
        
        xplotsize = 1.5
        if(ibmax ge 25) then xplotsize = 2.
        if(ibmax ge 50) then xplotsize = 2.5
        if(ibmax ge 100) then xplotsize = 3.
        if(ibmax ge 150) then xplotsize = 3.5
        
        yplotsize = 1.5

        mmax    = max(mass)
        mass(0) = mmax
        mmin    = min(mass)
        mleg    = fltarr(3)
        ileg    = intarr(3)
        if(mmin le 0.) then begin
           print,'fatal error mmin is nil'
           stop
        end
        ileg    = fltarr(3)
        ileg(2) = 0.6           ; symb size for mmin
        ileg(0) = 3.            ; symb size for mmax
        a       = (ileg(0)-ileg(2))/(alog10(mmax)-alog10(mmin))
        b       = (a*(alog10(mmax)+alog10(mmin))-(ileg(0)+ileg(2)))/(2.*a)
        mleg    = fltarr(3)
        ileg(1) = (ileg(0)+ileg(2))/2.
        for i = 0,2 do begin
           mleg(i) = 10^(ileg(i)/a + b)
        end
                
        
        file_eps = string(dir,"zoom_tree_",tree_number,"_",halo_number,"_",st_number,"_",stz_start,"_",stz_end,"_",method,".eps",format='(a,a,2(i8.8,a1),3(i3.3,a1),a3,a)')
        ;; xsize,ysize,!p.charsize,!x.charsize,!y.charsize,!p.charthick,feps,window,file
        init_plot,xplotsize,yplotsize,1.5,1.5,1.5,5.,feps,winnb,file_eps
        plot,[ibmin,ibmax],[stmin,stmax],xstyle=1,ystyle=1,/nodata,xmargin=[5,5],ymargin=[5,2]
        
        ;; plot squeleton of the merger tree
        for i = 1,n_elements(br_id) -1 do begin
           ison = son(i)
           if(ison gt 0) then begin
              oplot,[br_id(i),br_id(ison)],[st_id(i),st_id(ison)]
              if(f_lev(ison) gt 1 and f_lev(i) eq 1) then begin
                 ihost = host(ison)
                 if(ihost gt 0) then begin
                    oplot,[br_id(i),br_id(ihost)],[st_id(i),st_id(ihost)],linestyle=2
                 endif
              endif
           endif
        end
        
        for i = n_elements(br_id) -1,1,-1 do begin
           symb = a*(alog10(mass(i)) - b)
           icol = 0
           if(f_lev(i) gt 1) then begin
              symbols,30,0.25
              icol = 1
           endif else symbols,2,0.75
           if(f_frag(i) gt 0) then icol = 2
           plots,[br_id(i)],[st_id(i)],psym=8,symsize = symb + 1e-1
           plots,[br_id(i)],[st_id(i)],psym=8,symsize = symb - 1e-1,col=col(icol)
        end
        
        point_target,br_number,st_number,dst_leg,halo_number,col(2)
        if(feps eq 1) then begin
           device,/close
           set_plot,"x"
        endif
  
     endif
  end
  
end

;---------------------------------------------------------------------------------
pro find_targets,dir,tree_number,n_targ,halo_targnumb,st_targnumb,br_targnumb,st_targstart,st_targend

  read_where_haloes_are,dir,tree_numb,halo_numb,st_numb,branch_numb,branch_id,st_start,st_end

  n_targ = 0
  for itree = 0, n_elements(tree_numb)-1 do begin
     if(tree_numb(itree) eq tree_number) then n_targ = n_targ + 1
  end

  if(n_targ le 0) then return
  
  halo_targnumb = lonarr(n_targ)
  st_targnumb   = intarr(n_targ)
  br_targnumb   = intarr(n_targ)
  st_targend    = intarr(n_targ)
  st_targstart  = intarr(n_targ)
  
  i_targ = 0
  for itree = 0, n_elements(tree_numb)-1 do begin
     if(tree_numb(itree) eq tree_number) then begin
        halo_targnumb(i_targ) = halo_numb(itree)
        st_targnumb(i_targ)   = st_numb(itree)
        br_targnumb(i_targ)   = branch_numb(itree)
        st_targstart(i_targ)  = st_start(itree)
        st_targend(i_targ)    = st_end(itree)
        i_targ                = i_targ + 1
     endif
  end

end

;---------------------------------------------------------------------------------
pro point_target,xtarg,ytarg,dy,halo_targ,colid
     
  symbols,1,1.
  plots,[xtarg],[ytarg],psym=8,symsize = 4,col=colid,thick=2.
  plots,[xtarg],[ytarg],psym=1,symsize = 5,col=colid,thick=2.
  xyouts,xtarg,ytarg+dy,string(ytarg,halo_targ,format="(i3,1x,i8)"),alignment=0.5,col=colid

end

;---------------------------------------------------------------------------------
pro read_where_haloes_are,dir,tree_numb,halo_numb,st_numb,branch_numb,branch_id,st_start,st_end

  wherefile = string(dir,"where_haloes_are.dat",format='(a,a)')

  if(file_test(wherefile) ne 1) then begin
     print,"file ",wherefile," doesn't exist check var dir"
     stop
  endif
  if(file_test("wh_template.sav") eq 1) then begin
     restore,"wh_template.sav"
  endif else begin
     if(n_elements(template_wh) le 0) then template_wh = ascii_template(wherefile)
     save,template_wh,file="wh_template.sav"
  endelse
  data_wh = read_ascii(wherefile,template=template_wh)

  tree_numb   = data_wh.field1
  halo_numb   = data_wh.field2
  st_numb     = data_wh.field3
  branch_numb = data_wh.field4
  branch_id   = data_wh.field5
  st_end      = data_wh.field6
  st_start    = data_wh.field7

end

;---------------------------------------------------------------------------------
pro read_tree,mt_file,br_id,st_id,son,host,mass,f_lev,f_frag

  if(file_test(mt_file) ne 1) then begin
     print,"file ",mt_file," doesn't exist"
     stop
  endif
  if(file_test("mt_template.sav") eq 1) then begin
     restore,"mt_template.sav"
  endif else begin
     if(n_elements(template_mt) le 0) then template_mt = ascii_template(mt_file)
     save,template_mt,file="mt_template.sav"
  endelse
  data_mt = read_ascii(mt_file,template=template_mt)

  br_id   = data_mt.field1
  st_id   = data_mt.field2
  son     = data_mt.field3
  host    = data_mt.field4
  mass    = data_mt.field5
  f_lev   = data_mt.field6
  f_frag  = data_mt.field7
  
end

;---------------------------------------------------------------------------------
; xsize,ysize,!p.charsize,!x.charsize,!y.charsize,!p.charthick,feps,window,file
pro init_plot,xs,ys,csize,xcsize,ycsize,cth,feps,winnb,file_eps
  
  ;; auto parameter for output, X11 or eps
  ;; xs: image x size, ys; image y size
  ;; csize: plot charsize
  ;; xcsize: xaxis legend charsize, ysize: yaxis charsize
  ;; cth: character and lines thickness
  ;; feps: 1 to ouput eps file, 0 to display in X11
  ;; winnb: X11 window id
  ;; file_eps: name of the eps output file 
  

  if(feps ne 1) then begin
     !p.charsize  = csize
     !p.charthick = 1
     !p.thick     = 1
     !x.thick     = 1
     !y.thick     = 1
     !x.charsize  = xcsize
     !y.charsize  = ycsize
     if(winnb ge 0) then begin
        winid = winnb
        winnb = winnb + 1
     endif else begin
        winid = 0
     end
     set_plot,'x'
     window,winid,xsize = 300*xs,ysize=300*ys,retain=2
     device, decomposed = 0
  endif else begin
     !p.charsize  = csize
     !p.charthick = cth
     !p.thick     = cth
     !x.thick     = cth
     !y.thick     = cth
     !x.charsize  = xcsize
     !y.charsize  = ycsize
     set_plot, 'ps'
     device,filename = file_eps, xsize = 10.*xs, ysize = 10.*ys,/color,bits_per_pixel=8,/portrait, /encapsulated,scale_factor=0.5
     print, 'exporting: ',file_eps
  endelse
  
end
