;-----------------------------------
pro plot_all_list,dir,method,feps

  winnb = 0
  
  ; get_merger_tree_list
  get_merger_tree_list,dir,tree_list
  
  for ilist = 0,n_elements(tree_list) -1 do begin
     plot_full_merger_tree,dir,tree_list(ilist),method,feps,winnb
  end
  
  
end

;-----------------------------------
pro plot_full_merger_tree,dir,tree_number,method,feps,winnb


  loadct, 12,ncolors=16
  file_tree = string(dir,"full_tree_",tree_number,"_",method,".dat",format='(a,a,i8.8,a1,a3,a)')
  read_merger_tree,file_tree,ib,st,is,ih,m,fs,ff
  
  ibmax = max(ib)
  ib(0) = ibmax
  ibmin = - ibmax/50
  ibmax = ibmax - (ibmin - min(ib)) 
 
  stmax = max(st)
  st(0) = stmax
  stmin = min(st)
  stmin = stmin - stmin mod 10

  stup   = (stmax-stmin)/20.
;  stmax  = stmax + stup
  
  
  col    = [6,7,12]
  
  file_eps = string(dir,"full_merger_tree_",tree_number,"_",method,".eps",format='(a,a,i8.8,a1,a3,a)')
  xpsize   = 2.
  if(ibmax ge 100) then xpsize = 2.5
  if(ibmax ge 250) then xpsize = 3.
  if(ibmax ge 500) then xpsize = 4.
  if(ibmax ge 750) then xpsize = 4.5
  if(ibmax ge 1000) then xpsize = 5.
  if(ibmax ge 1250) then xpsize = 5.5
  if(ibmax ge 1500) then xpsize = 6.
  ypsize   = 2.
  if(stmax ge 50) then ypsize = 2.5
  if(stmax ge 75) then ypsize = 3.
  if(stmax ge 100) then ypsize = 3.5 
  
  print,"xpsize,ypsize",xpsize,ypsize
  ;; xsize,ysize,!p.charsize,!x.charsize,!y.charsize,!p.charthick,feps,window,file
  init_plot,xpsize,ypsize,1.5,1.,1.,7,feps,winnb,file_eps
  !p.multi = [0,1,1]

  plot,[ibmin,ibmax],[stmin,stmax],xstyle=1,ystyle = 1,/nodata,ymargin=[4,4],xmargin=[5,2]
;;   axis,xaxis=0,xrange = [ibmin,ibmax],xstyle = 1
;;   axis,yaxis=0,yrange = [stmin,stmax],ystyle = 1
;;   axis,yaxis=1,yrange = [stmin,stmax],ystyle = 1
  for i = 1L,n_elements(ib) -1 do begin
     ison = is(i)
     if(ison gt 0) then begin
        icol = 0 ; default  halos
        if(fs(i) gt 1) then icol = 1 ; subhalo
        if(ff(i) gt 0) then icol = 2 ; fragment
        oplot,[ib(i),ib(ison)],[st(i),st(ison)],col=col(icol)
;;         if(ib(ison) eq ib(i)) then begin
;;            if(fs(ison) gt 1) then begin
;;               oplot,[ib(i),ib(ison)],[st(i),st(ison)],linestyle=2
;;            endif else begin
;;               oplot,[ib(i),ib(ison)],[st(i),st(ison)]
;;            endelse
;;         endif else begin
;;            oplot,[ib(i),ib(ison)],[st(i),st(ison)]
;;         endelse
     endif else begin
        ihost = ih(i)
   
        if(st(i) eq st(1) and ihost gt 0) then begin
           plots,[ib(ihost),ib(ihost)],[st(ihost),st(ihost)+stup]
           plots,[ib(i),ib(ihost)],[st(i),st(ihost)+stup]
        end
     endelse
     if(ff(i) eq 1) then oplot,[ib(i)],[st(i)],psym=7 
  end
  get_location,dir,tree_number,nloc,idhalo,sthalo,ibhalo,ibloc,ststartloc,stendloc
  
  if(nloc gt 0) then begin
     for iloc = 0,nloc-1 do begin
        plot_target,ibhalo(iloc),sthalo(iloc),stup,idhalo(iloc),col(2)
     end
  end
  if(feps eq 1) then begin
     device,/close
     set_plot,"x"
  end
  if(nloc gt 0) then begin
     for iloc = 0,nloc- 1 do begin
        print,"plot zoom",sthalo(iloc),idhalo(iloc),ststartloc(iloc),stendloc(iloc),format="(a,i3,1x,i8,2(1x,i3))"
        plot_zoom_merger_tree,dir,tree_number,idhalo(iloc),sthalo(iloc),ibloc(iloc),ststartloc(iloc),stendloc(iloc),method,feps,winnb
     end
  endif

end

;-----------------------------------
pro plot_zoom_merger_tree,dir,tree_number,halo_number,step_number,ib_loc,ststart,stend,method,feps,winnb

  file_tree = string(dir,"zoom_tree_",tree_number,"_",halo_number,"_",step_number,"_",ststart,"_",stend,"_",method,".dat",format='(a,a,2(i8.8,a1),3(i3.3,a1),a3,a)')
  
  read_merger_tree,file_tree,ib,st,is,ih,m,fs,ff
  
  if(n_elements(m) le 2) then begin
     print,'only one halo in merger treen, returning'
     return
  end
 
  ibmax = max(ib)
  ib(0) = ibmax
  ibmin = 0
  ibmax = ibmax + (min(ib)-ibmin)
  
  stmax = max(st)
  st(0) = stmax
  stmin = min(st)
  stup  = (stmax-stmin)/20.
  stleg = stmin - (stmax-stmin)/10.
  xleg  = fltarr(3,3)
  for ileg = 0,2 do begin 
     xleg(ileg,2) = ibmin + (ibmax-ibmin)*((2.*ileg+0.5))/6.
     xleg(ileg,1) = xleg(ileg,2) - (ibmax-ibmin)/((ileg+1)*18.)
     xleg(ileg,0) = xleg(ileg,2) - (ibmax-ibmin)/((ileg+1)*9.)
  end
  
                             
; mass legend
  mmax        = max(m)
  m(0)        = mmax
  mmin        = min(m)
    
  ilegsize    = fltarr(3)
  ilegsize(2) = 0.2                 ; symb size for mmin
  ilegsize(0) = 1.                  ; symb size for mmax
  a           = (ilegsize(0)-ilegsize(2))/(alog10(mmax)-alog10(mmin))
  b           = (a*(alog10(mmax)+alog10(mmin))-(ilegsize(0)+ilegsize(2)))/(2.*a)
  mleg        = fltarr(3)
  mleglog     = fltarr(3)
  ilegsize(1) = (ilegsize(0)+ilegsize(2))/2.
  for ileg = 0,2 do begin
     mleg(ileg)    = 10^(ilegsize(ileg)/a + b)
     mleglog(ileg) = alog10(mleg(ileg)) - alog10(mleg(ileg)) mod 1
     mleg(ileg)    = 10.^(alog10(mleg(ileg)) mod 1)
  end
  symb    = fltarr(n_elements(m))


  col = [6,12,7,12]
  file_eps = string(dir,"zoom_tree_",tree_number,"_",halo_number,"_",step_number,"_",ststart,"_",stend,"_",method,".eps",format='(a,a,2(i8.8,a1),3(i3.3,a1),a3,a)')
  ;; xsize,ysize,!p.charsize,!x.charsize,!y.charsize,!p.charthick,feps,window,file
  init_plot,2.,2.,1.5,1.,1.,7,feps,winnb,file_eps
  !p.multi = [0,1,1]
  ; mass legend
  symb    = fltarr(n_elements(m))
  
  plot,[ibmin,ibmax],[stmin,stmax],xstyle=1,ystyle = 1,/nodata,xmargin=[4,4],ymargin=[6,4]
  ; squeleton of the merger tree
  for i = 1L,n_elements(ib) -1 do begin
     ison = is(i)
     if(ison gt 0) then begin
        oplot,[ib(i),ib(ison)],[st(i),st(ison)]
        if(fs(ison) gt 1 and fs(i) le 1) then begin
           ihost = ih(ison)
           if(ihost gt 0) then begin
              oplot,[ib(i),ib(ihost)],[st(i),st(ihost)],linestyle=2
           endif
        endif
     endif
     symb(i) = a*(alog10(m(i)) - b) 
  end  
  for i = 1L,n_elements(ib) -1 do begin
     icol = 0
     if(fs(i) gt 1) then begin
        symbols,30,1.
        icol = icol + 2
     end else symbols,2,3.
     if(ff(i) gt 0) then icol = icol + 1
     plots,[ib(i)],[st(i)],psym=8,symsize=symb(i)+0.05
     plots,[ib(i)],[st(i)],psym=8,symsize=symb(i)-0.05,col=col(icol)
  end
  if(ib_loc gt 0) then plot_target,ib_loc,step_number,stup,halo_number,col(1)
 
  ; legend
  symbols,2,3.
  for ileg = 0,2 do begin
     plots,[xleg(ileg,0)],[stleg],psym=8,symsize=ilegsize(ileg)+0.05
     plots,[xleg(ileg,0)],[stleg],psym=8,symsize=ilegsize(ileg)-0.05,col=col(0)
  end
  
  symbols,30,1.
  for ileg = 0,2 do begin
     plots,[xleg(ileg,1)],[stleg],psym=8,symsize=ilegsize(ileg)+0.05
     plots,[xleg(ileg,1)],[stleg],psym=8,symsize=ilegsize(ileg)-0.05,col=col(2)
  end  

  for ileg = 0,2 do begin
     xyouts,xleg(ileg,2),stleg,string(mleg(ileg)," 10!u",mleglog(ileg),"!n M!D!Mn!X!N",format='(F4.2,a,i2.2,a)')
  end

  if(feps eq 1) then begin
     device,/close
     set_plot,"x"
  end
  
end

;-----------------------------------
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
     device,filename = file_eps, xsize = 10.*xs, ysize = 10.*ys,/color,bits_per_pixel=8,/portrait, /encapsulated;,scale_factor=0.5
     print, 'exporting: ',file_eps
  endelse
  
end

;-----------------------------------
pro get_location,dir,tree_number,nloc,idhalo,sthalo,ibhalo,ibloc,ststloc,stendloc
  
  read_where_halo_are,dir,t_numb,h_numb,s_numb,b_numb,b_id,s_end,s_start
  
  nloc = 0
  for itree = 0,n_elements(t_numb)-1 do begin
     if(t_numb(itree) eq tree_number) then nloc = nloc + 1
  end
  if(nloc le 0) then return
  

  sthalo   = intarr(nloc)
  idhalo   = fltarr(nloc)
  ibhalo   = fltarr(nloc)
  ibloc    = fltarr(nloc)
  ststloc  = intarr(nloc)
  stendloc = intarr(nloc)
  nloc = 0
  for itree = 0,n_elements(t_numb)-1 do begin
     if(t_numb(itree) eq tree_number) then begin
        sthalo(nloc)   = s_numb(itree)
        idhalo(nloc)   = h_numb(itree)
        ibhalo(nloc)   = b_numb(itree)
        ibloc(nloc)    = b_id(itree)
        ststloc(nloc)  = s_start(itree)
        stendloc(nloc) = s_end(itree)
        nloc           = nloc + 1
     endif
  end
  
end

;-----------------------------------
pro get_merger_tree_list,dir,tree_list

  read_where_halo_are,dir,t_numb,h_numb,s_numb,b_numb,b_id,s_end,s_start
  
  ntree = 1
  if(n_elements(t_numb) gt 1) then begin
     for i = 1,n_elements(t_numb) - 1 do begin
        ffound  = 0
        tt_numb = t_numb(i) 
        for j = 0,i-1 do begin
           if(t_numb(j) eq tt_numb) then ffound = 1
        end
        if(ffound eq 0) then ntree = ntree + 1
     end
  endif
  
  tree_list        = intarr(ntree)
  ntree            = 0
  tree_list(ntree) = t_numb(0)
  if(n_elements(t_numb) gt 1) then begin
     for i = 1,n_elements(t_numb) - 1 do begin
        ffound  = 0
        tt_numb = t_numb(i) 
        for j = 0,i-1 do begin
           if(t_numb(j) eq tt_numb) then ffound = 1
        end
        if(ffound eq 0) then begin
           ntree = ntree + 1
           tree_list(ntree) = t_numb(i)
        endif
     end
  endif 

end

;-----------------------------------
pro plot_target,x_targ,y_targ,y_up,halo,colred
 
  symbols,1,1
  plots,[x_targ],[y_targ],psym=1,symsize = 5,col=colred,thick=2*!p.thick
  plots,[x_targ],[y_targ],psym=8,symsize = 3,col=colred,thick=4*!p.thick
  xyouts,x_targ,y_targ+y_up,string(y_targ,halo,format='(i3,i8)'),col=colred,alignment=0.5

end
;-----------------------------------
pro read_merger_tree,file_tree,ib,st,is,ih,m,fs,ff
  
  ifexist = file_test('mt_template.sav')
  if(ifexist eq 1) then begin
     restore,'mt_template.sav'
  endif else begin
     if(n_elements(template_mt) eq 0) then template_mt = ascii_template(file_tree)
     save,template_mt,file="mt_template.sav"
  endelse
  print, 'now reading ',file_tree
  data_mt = read_ascii(file_tree, template=template_mt)
  ib      = data_mt.field1
  st      = data_mt.field2
  is      = data_mt.field3
  ih      = data_mt.field4
  m       = data_mt.field5
  fs      = data_mt.field6
  ff      = data_mt.field7

end

;-----------------------------------
pro read_where_halo_are,dir,t_numb,h_numb,s_numb,b_numb,b_id,s_end,s_start

  wherefile = string(dir,"where_haloes_are.dat",format='(a,a)')
  ifexist = file_test('wh_template.sav')
  if(ifexist eq 1) then begin
     restore,'wh_template.sav'
  endif else begin
     if(n_elements(template_wh) eq 0) then template_wh = ascii_template(wherefile)
     save,template_wh,file="wh_template.sav"
  endelse 

  data_wh = read_ascii(wherefile, template=template_wh)
  t_numb  = data_wh.field1
  h_numb  = data_wh.field2
  s_numb  = data_wh.field3
  b_numb  = data_wh.field4
  b_id    = data_wh.field5
  s_end   = data_wh.field6
  s_start = data_wh.field7
  
end

