function main_branch,haloID,ts,t,p,sim

  ;; select tree rooted on haloID 
  ihalo = where(t.haloid eq haloID,n)
  if n eq 0 then stop,'HaloID not found in main_branch...'
  lastID = t(ihalo).lastProgenitorID
  sel    = where(t.haloid ge haloID and t.haloid le lastID)  ;; all progs of selected halo
  tree   = t(sel)   
  props  = p(sel) 
  
  ;; sort tree with haloid 
  ind   = sort(tree.haloid) 
  tree  = tree(ind)
  props = props(ind)

  ;; get branches
  leaves = where(tree.haloid eq tree.lastProgenitorID)       ;; end of branches 
  lastID = tree(leaves(0)).haloid  ;; end of first branch (i.e. main branch) 
  branch = where(tree.haloid ge haloID and tree.haloid le lastID,nbranch)

  ;; transform branch back to global index ...   
  branch = sel(ind(branch))
  
  return,branch

end


function read_timesteps_props,file

  nsteps = (r = 0L) 
  openr,11,file
  readu,11,r,nsteps,r
  
  ts = replicate({nh:0L,aexp:0.0d0,age_univ:0.0d0},nsteps)
  nh = lonarr(nsteps)
  readu,11,r,nh,r
  ts(*).nh = nh
  x = fltarr(nsteps)
  readu,11,r,x,r
  ts(*).aexp = x
  readu,11,r,x,r
  ts(*).age_univ = x
  close,11

  return,ts 

end


function read_tree_struct,file
  
  nsteps = (nIDs = (nIndex = (r = 0L)))
  openr,11,file
  readu,11,r,nsteps,nIDs,nIndex,r
  nh = lonarr(nsteps)
  readu,11,r,nh,r
  nhtot = total(nh)
  zero = long64(0)
  IDs = replicate({BushID:zero,TreeID:zero,HaloID:zero,haloNum:zero,$
                   haloTimestep:zero,FirstProgenitorID:zero,NextProgenitorID:zero,$
                   DescendentID:zero,LastProgenitorID:zero,HostHaloID:zero,$
                   HostSubID:zero,NextSubID:zero},nhtot)
  ihalo = 0L
  for ts = 0,nsteps-1 do begin 
     if nh(ts) gt 0 then begin 
        blah = lon64arr(nids,nh(ts))
        readu,11,r,blah,r
        for i = 0L,nIDs-1L do begin 
           ids(ihalo:ihalo+nh(ts)-1L).(i) = reform(blah(i,*))
        endfor
        ihalo = ihalo + nh(ts)
        ;; skip indexes
        blah = lonarr(nindex,nh(ts))
        readu,11,r,blah,r
     endif
  endfor
  close,11

  return,ids

end


function read_halo_props,file
  
  nsteps = (nprops = (r = 0L))
  openr,11,file
  readu,11,r,nsteps,nprops,r
  nh = lonarr(nsteps)
  readu,11,r,nh,r
  nhtot = total(nh)
  p = replicate({x:0.0,y:0.0,z:0.0,vx:0.0,vy:0.0,vz:0.0,mtot:0.0,r:0.0,spin:0.0,$
                 rvir:0.0,mvir:0.0,tvir:0.0,cvel:0.0,dmacc:0.0,frag:0.0,lx:0.0,$
                 ly:0.0,lz:0.0,ep:0.0,ek:0.0,et:0.0},nhtot)
  ihalo = 0L
  for ts = 0,nsteps-1 do begin 
     if nh(ts) gt 0 then begin 
        blah = fltarr(nprops,nh(ts))
        readu,11,r,blah,r
        for i = 0L,nprops-1L do begin 
           p(ihalo:ihalo+nh(ts)-1L).(i) = reform(blah(i,*))
        endfor
        ihalo = ihalo + nh(ts)
     endif
  endfor
  close,11

  return,p

end


function define_zoom,hid,ts,t,p,sim,file=file
  
  ;; INPUTS
  ;; 
  ;; - hid : the haloID of the halo to zoom on
  ;; 
  ;; - ts,t,p,sim : the tree and simulation.
  ;; 
  ;; KEYWORDS
  ;; 
  ;; - file : if set, output info to a file.

  i        = where(t.haloid eq hid)
  hbushid  = t(i).bushid
  halonum  = t(i).halonum 
  halostep = t(i).halotimestep
  xh       = p(i).x & yh = p(i).y & zh = p(i).z            ;; position of selected halo
  coords   = {x:xh/sim.lbox+0.5,y:yh/sim.lbox+0.5,z:zh/sim.lbox+0.5} ;; coords of selected halo in box units. 
  
  ;; select all halos in bush
  sel  = where(t.bushid eq hbushid)
  psel = p(sel) & tsel = t(sel)
  
  ;; recenter on selected halo and correct for periodic boundaries
  psel.x = psel.x - xh & psel.y = psel.y - yh & psel.z = psel.z - zh
  ii = where(psel.x gt 0.5*sim.lbox,ni)  & if ni ne 0 then psel(ii).x = psel(ii).x - sim.lbox
  ii = where(psel.x lt -0.5*sim.lbox,ni) & if ni ne 0 then psel(ii).x = psel(ii).x + sim.lbox
  ii = where(psel.y gt 0.5*sim.lbox,ni)  & if ni ne 0 then psel(ii).y = psel(ii).y - sim.lbox
  ii = where(psel.y lt -0.5*sim.lbox,ni) & if ni ne 0 then psel(ii).y = psel(ii).y + sim.lbox
  ii = where(psel.z gt 0.5*sim.lbox,ni)  & if ni ne 0 then psel(ii).z = psel(ii).z - sim.lbox
  ii = where(psel.z lt -0.5*sim.lbox,ni) & if ni ne 0 then psel(ii).z = psel(ii).z + sim.lbox

  ;; define range (re-center on center, not on final halo)
  xmin = min(psel.x) & xmax = max(psel.x) & xc = 0.5 * (xmin + xmax) & dx = xmax - xmin
  ymin = min(psel.y) & ymax = max(psel.y) & yc = 0.5 * (ymin + ymax) & dy = ymax - ymin
  zmin = min(psel.z) & zmax = max(psel.z) & zc = 0.5 * (zmin + zmax) & dz = zmax - zmin
  
  ;; back to original coords
  xc = xc + xh & yc = yc + yh & zc = zc + zh
  ;; and box units (-> from 0 to 1)
  xc = xc / sim.lbox + 0.5 & yc = yc / sim.lbox + 0.5 & zc = zc / sim.lbox + 0.5
  dx = dx / sim.lbox & dy = dy / sim.lbox & dz = dz / sim.lbox
  
  zoom = {x:xc,y:yc,z:zc, $
          x128:round(xc*128),y128:round(yc*128),z128:round(zc*128), $
          x256:round(xc*256),y256:round(yc*256),z256:round(zc*256), $
          x512:round(xc*512),y512:round(yc*512),z512:round(zc*512), $
          x1024:round(xc*1024),y1024:round(yc*1024),z1024:round(zc*1024), $
          x2048:round(xc*2048),y2048:round(yc*2048),z2048:round(zc*2048), $
          dx:dx,dy:dy,dz:dz}

  if keyword_set(file) then begin 
     openw,11,file
     printf,11,'# informations to zoom on halo '+strtrim(halonum,2)+' at (tree) timestep '+strtrim(halostep,2)
     printf,11,'#'
     printf,11,'# center in box units (from 0 to 1)'
     printf,11,zoom.x,zoom.y,zoom.z
     printf,11,'# center in pixel units: '
     printf,11,'@128: ',zoom.x128,zoom.y128,zoom.z128
     printf,11,'@256: ',zoom.x256,zoom.y256,zoom.z256
     printf,11,'@512: ',zoom.x512,zoom.y512,zoom.z512
     printf,11,'@1024: ',zoom.x1024,zoom.y1024,zoom.z1024
     printf,11,'@2048: ',zoom.x2048,zoom.y2048,zoom.z2048
     printf,11,'# dx, dy, dz :'
     printf,11,'# BOX UNITS: ',zoom.dx,zoom.dy,zoom.dz
     printf,11,'@128:',round(zoom.dx*128),round(zoom.dy*128),round(zoom.dz*128)
     printf,11,'@256:',round(zoom.dx*256),round(zoom.dy*256),round(zoom.dz*256)
     printf,11,'@512:',round(zoom.dx*512),round(zoom.dy*512),round(zoom.dz*512)
     printf,11,'@1024:',round(zoom.dx*1024),round(zoom.dy*1024),round(zoom.dz*1024)
     printf,11,'@2048:',round(zoom.dx*2048),round(zoom.dy*2048),round(zoom.dz*2048)
     close,11
  endif

  return, zoom

end


pro read_tree,dir,snapnum,ts,t,p,sim
  
  ;; INPUTS:
  ;; 
  ;; - dir : where files outputed by TreeMaker are. 
  ;; - snapnum : the last timestep used to generate trees (appears in
  ;;   filename). 
  ;; 

  file = dir+'/tstep_file_'+string(snapnum,format="(i3.3)")+'.001' ;; '../ref/tstep_file_018.001'
  ts = read_timesteps_props(file)
  file = dir+'/tree_file_'+string(snapnum,format="(i3.3)")+'.001'  ;; '../ref/tree_file_018.001'
  t = read_tree_struct(file)
  file = dir+'/props_'+string(snapnum,format="(i3.3)")+'.001'      ;;'../ref/props_018.001'
  p = read_halo_props(file)

  ;; define simulation properties
  sim = {Lbox:20.0,h:0.702,Om:0.277,Ol:0.723}
  
  ;; convert postitions to comoving coords in Mpc/h
  st = t.halotimestep - 1
  aexp = ts(st).aexp
  p.x = p.x / aexp * sim.h
  p.y = p.y / aexp * sim.h
  p.z = p.z / aexp * sim.h

end


pro plot_all_trees,dir,snapnum, $
                   psfile=psfile, half_plot_size=half_plot_size, $
                   minmass=minmass,maxmass=maxmass

  ;; INPUTS:
  ;; 
  ;; - dir : where files outputed by TreeMaker are. 
  ;; - snapnum : the last timestep used to generate trees (appears in
  ;;   filename). 
  ;; 

  if keyword_set(psfile) then begin 
     set_plot,'ps'
     device,filename=psfile,/color,bits_per_pixel=8,xsize=20,ysize=20,yoffset=3,xoffset=0.5
  endif
  
  ;; read tree files and define simulation properties.
  read_tree,dir,snapnum,ts,t,p,sim
  
  ;; SELECTIONS ---- 
  
  ;; massive halos at z = 3
  z     = 3.
  aexp  = 1./(1+z) 
  da    = abs(ts.aexp - aexp)
  tsnum = (where(da eq min(da)))(0) + 1 ;; timestep (relative to tree, NOT output nb)
  tsnum = 16
  z3 = where(t.halotimestep eq tsnum)
  if not keyword_set(minmass) then minmass = 5.
  if not keyword_set(maxmass) then maxmass = 1000000.
  hlist = where(p(z3).mtot gt minmass and p(z3).mtot lt maxmass)
  haloidlist = t(z3(hlist)).haloid

  for ihalo = 0,n_elements(haloidlist)-1 do begin 
     haloid = haloidlist(ihalo)
     plot_halo_tree,haloid,ts,t,p,sim,scale='mass',half_plot_size=half_plot_size,/highlight_bush
  endfor

  ;; list of all treeIDs 
  ;; list = t.treeID 
  ;; list = list(uniq(sort(list)))
  ;; for itree = 0L,n_elements(list)-1 do begin 
  ;;    tt = where(t.treeid eq list(itree),np)
  ;;    if np gt 200 then begin 
  ;;       print,itree
  ;;       haloid = min(t(tt).haloid)
  ;;       plot_halo_tree,haloid,ts,t,p
  ;;    endif
  ;; endfor

  if keyword_set(psfile) then begin 
     device,/close
     set_plot,'x'
  endif

end
 

pro plot_halo_tree,fid,ts,t,p,sim,scale=scale,half_plot_size=half_plot_size, $
                   show_subhalos=show_subhalos, highlight_bush=highlight_bush
  
  ;; INPUTS
  ;;
  ;; - fid : HaloID of halo for which we should plot tree. 
  ;; 
  ;; - ts  : timestep information (from read_timesteps_props)
  ;; 
  ;; - t   : tree information (from read_tree_struct)
  ;;
  ;; - p   : halo properties (from read_halo_props)
  ;;
  ;; - sim : simulation info (hardcoded)
  ;; 
  ;; 
  ;; KEYWORDS 
  ;;
  ;; - highlight_bush : if set, will highlight all halos which belong
  ;;   to the same bush as selected halo.
  ;; 

  ;; fid is haloid of root of tree -> get index of that halo in list
  ;; and extract some useful props
  i        = where(t.haloid eq fid)
  hbushid  = t(i).bushid
  halonum  = t(i).halonum 
  halostep = t(i).halotimestep
  xh       = p(i).x & yh = p(i).y & zh = p(i).z            ;; position of selected halo
  coords   = {x:xh/sim.lbox+0.5,y:yh/sim.lbox+0.5,z:zh/sim.lbox+0.5} ;; coords of selected halo in box units. 
  lid      = t(i).lastprogenitorid                         ;; ID of last progenitor
  sel      = where(t.haloid ge fid and t.haloid le lid,np) ;; indexes of progenitors.

  ;; substract offset to IDs of selected tree, and sort them ascending
  ;; so the can be used for local indexing
  psel = p(sel) 
  tsel = t(sel) 
  ind  = sort(tsel.haloid)
  psel = psel(ind)
  tsel = tsel(ind) 
  tsel.haloid = tsel.haloid - fid
  tsel.descendentid = tsel.descendentid - fid
  
  ;; plots
  window,retain=2,xsize=1000,ysize=1000
  device,decomposed=0
  loadct,0,/silent
  !p.multi=[0,2,2]

  ;; recenter on selected halo and correct for periodic boundaries
  psel.x = psel.x - xh & psel.y = psel.y - yh & psel.z = psel.z - zh
  ii = where(psel.x gt 0.5*sim.lbox,ni)  & if ni ne 0 then psel(ii).x = psel(ii).x - sim.lbox
  ii = where(psel.x lt -0.5*sim.lbox,ni) & if ni ne 0 then psel(ii).x = psel(ii).x + sim.lbox
  ii = where(psel.y gt 0.5*sim.lbox,ni)  & if ni ne 0 then psel(ii).y = psel(ii).y - sim.lbox
  ii = where(psel.y lt -0.5*sim.lbox,ni) & if ni ne 0 then psel(ii).y = psel(ii).y + sim.lbox
  ii = where(psel.z gt 0.5*sim.lbox,ni)  & if ni ne 0 then psel(ii).z = psel(ii).z - sim.lbox
  ii = where(psel.z lt -0.5*sim.lbox,ni) & if ni ne 0 then psel(ii).z = psel(ii).z + sim.lbox

  ;; define plot range (re-center on center, not on final halo)
  xmin = min(psel.x) & xmax = max(psel.x) & xc = 0.5 * (xmin + xmax)
  ymin = min(psel.y) & ymax = max(psel.y) & yc = 0.5 * (ymin + ymax)
  zmin = min(psel.z) & zmax = max(psel.z) & zc = 0.5 * (zmin + zmax)
  if keyword_set(half_plot_size) then dd = half_plot_size else dd = 1.0d0 
  xmin = (ymin = (zmin = -dd)) & xmax = (ymax = (zmax = dd))
  psel.x = psel.x - xc
  psel.y = psel.y - yc
  psel.z = psel.z - zc

  ;; re-center other halos too ... and correct periodic boundaries
  x = p.x - xc - xh 
  ii = where(x gt 0.5*sim.lbox,ni)  & if ni ne 0 then x(ii) = x(ii) - sim.lbox
  ii = where(x lt -0.5*sim.lbox,ni) & if ni ne 0 then x(ii) = x(ii) + sim.lbox
  y = p.y - yc - yh
  ii = where(y gt 0.5*sim.lbox,ni)  & if ni ne 0 then y(ii) = y(ii) - sim.lbox
  ii = where(y lt -0.5*sim.lbox,ni) & if ni ne 0 then y(ii) = y(ii) + sim.lbox
  z = p.z - zc - zh
  ii = where(z gt 0.5*sim.lbox,ni)  & if ni ne 0 then z(ii) = z(ii) - sim.lbox
  ii = where(z lt -0.5*sim.lbox,ni) & if ni ne 0 then z(ii) = z(ii) + sim.lbox

  ;; select other halos within plot range
  others = where(x lt dd and x gt -dd and $
                 y lt dd and y gt -dd and $
                 z lt dd and z gt -dd,nothers)
  if nothers gt 0 then begin 
     osc = alog10(p(others).mtot)
     osc = (osc - min(osc)) / (max(osc) - min(osc)) * 4.
  endif else begin 
     osc = [0]
  endelse

  ;; symbol sizes 
  if keyword_set(scale) then begin 
     if scale eq 'mass' then begin 
        sc = alog10(psel.mtot)
        sc = (sc - min(sc)) / (max(sc)-min(sc)) * 4.
     endif
  endif else begin 
     sc = psel.mtot * 0.0 + 1.0
  endelse
  ;; symbol colors
  ;; selected tree halos
  couleur = bytscl(tsel.halotimestep,top=254)
  ;; other halos
  ocouleur = bytscl(t.halotimestep,top=160) + 40   ;; fill-color is timestep
  oct      = ocouleur * 0    ;; default color-table for other halos (grey scale)
  if keyword_set(highlight_bush) then begin 
     ;; change circle color of bush halos
     ii = where(t.bushid eq hbushid)
     oct(ii) = 38 ;; will turn color 180 into red... 
  endif

  ;; x vs y
  loadct,0,/silent
  plot,psel.x,psel.y,psym=3,xtitle='x (cMpc)',ytitle='y (cMpc)',xr=[xmin,xmax],yr=[ymin,ymax],/nodata,title='timestep:'+strtrim(halostep,2)+', hno:'+strtrim(halonum,2)
  if nothers gt 0 then begin 
     for i = 0L,n_elements(others)-1 do begin 
        if keyword_set(highlight_bush) then loadct,0,/silent
        jeje_symbols,2,1.5
        plots,[x(others(i))],[y(others(i))],psym=8,col=ocouleur(others(i)),symsize=osc(i)
        if keyword_set(highlight_bush) then loadct,oct(others(i)),/silent
        jeje_symbols,1,1.5
        plots,[x(others(i))],[y(others(i))],psym=8,col=180,symsize=osc(i)
     end
  endif
  for i = 0,n_elements(sel)-1 do begin 
     if tsel(i).descendentid ge 0 then begin
        i2 = tsel(i).descendentid
        oplot,[psel(i).x,psel(i2).x],[psel(i).y,psel(i2).y],col=80
     endif
  endfor
  for i = 0,n_elements(sel)-1 do begin 
     loadct,39,/silent
     jeje_symbols,2,1.5
     plots,[psel(i).x],[psel(i).y],color=couleur(i),psym=8,symsize=sc(i)
     loadct,0,/silent
     jeje_symbols,1,1.5
     plots,[psel(i).x],[psel(i).y],color=80,psym=8,symsize=sc(i)
  endfor

  title='(x,y,z) = ('+strtrim(coords.x,2)+','+strtrim(coords.y,2)+','+strtrim(coords.z,2)+')  [ box units ]' + $
        '!C !C '+$
        ' @ 2048 = ('+strtrim(fix(coords.x*2048),2)+','+strtrim(fix(coords.y*2048),2)+','+strtrim(fix(coords.z*2048),2)+')'+$
        '   ,   '+$
        ' @ 1024 = ('+strtrim(fix(coords.x*1024),2)+','+strtrim(fix(coords.y*1024),2)+','+strtrim(fix(coords.z*1024),2)+')'+$
        '!C !C '+$
        ' @ 512  = ('+strtrim(fix(coords.x*512),2)+','+strtrim(fix(coords.y*512),2)+','+strtrim(fix(coords.z*512),2)+')'+$
        '   ,   '+$
        ' @ 256  = ('+strtrim(fix(coords.x*256),2)+','+strtrim(fix(coords.y*256),2)+','+strtrim(fix(coords.z*256),2)+')'
  xyouts,0,1.6*dd,title

  ;; x vs z
  plot,psel.x,psel.z,psym=3,xtitle='x (cMpc)',ytitle='z (cMpc)',xr=[xmin,xmax],yr=[zmin,zmax],/nodata
  if nothers gt 0 then begin 
     for i = 0L,n_elements(others)-1 do begin 
        if keyword_set(highlight_bush) then loadct,0,/silent
        jeje_symbols,2,1.5
        plots,[x(others(i))],[z(others(i))],psym=8,col=ocouleur(others(i)),symsize=osc(i)
        if keyword_set(highlight_bush) then loadct,oct(others(i)),/silent
        jeje_symbols,1,1.5
        plots,[x(others(i))],[z(others(i))],psym=8,col=180,symsize=osc(i)
     end
  endif
  for i = 0,n_elements(sel)-1 do begin 
     if tsel(i).descendentid ge 0 then begin
        i2 = tsel(i).descendentid
        oplot,[psel(i).x,psel(i2).x],[psel(i).z,psel(i2).z],col=80
     endif
  endfor
  for i = 0,n_elements(sel)-1 do begin 
     loadct,39,/silent
     jeje_symbols,2,1.5
     plots,[psel(i).x],[psel(i).z],color=couleur(i),psym=8,symsize=sc(i)
     loadct,0,/silent
     jeje_symbols,1,1.5
     plots,[psel(i).x],[psel(i).z],color=80,psym=8,symsize=sc(i)
  endfor

  ;; y vs z
  plot,psel.y,psel.z,psym=3,xtitle='y (cMpc)',ytitle='z (cMpc)',xr=[ymin,ymax],yr=[zmin,zmax],/nodata
  if nothers gt 0 then begin 
     for i = 0L,n_elements(others)-1 do begin 
        if keyword_set(highlight_bush) then loadct,0,/silent
        jeje_symbols,2,1.5
        plots,[y(others(i))],[z(others(i))],psym=8,col=ocouleur(others(i)),symsize=osc(i)
        if keyword_set(highlight_bush) then loadct,oct(others(i)),/silent
        jeje_symbols,1,1.5
        plots,[y(others(i))],[z(others(i))],psym=8,col=180,symsize=osc(i)
     end
  endif
  for i = 0,n_elements(sel)-1 do begin 
     if tsel(i).descendentid ge 0 then begin
        i2 = tsel(i).descendentid
        oplot,[psel(i).y,psel(i2).y],[psel(i).z,psel(i2).z],col=80
     endif
  endfor
  for i = 0,n_elements(sel)-1 do begin 
     loadct,39,/silent
     jeje_symbols,2,1.5
     plots,[psel(i).y],[psel(i).z],color=couleur(i),psym=8,symsize=sc(i)
     loadct,0,/silent
     jeje_symbols,1,1.5
     plots,[psel(i).y],[psel(i).z],color=80,psym=8,symsize=sc(i)
  endfor

  ;; plot mass history
  plot_io,tsel.halotimestep,psel.mtot,psym=3,ytitle='M!dh!n [10!u11!n M!d!9n!n!3]',xtitle='ts',title='HaloID = '+strtrim(fid,2)
  loadct,39,/silent
  jeje_symbols,2,1.5
  plots,tsel.halotimestep,psel.mtot,psym=8,color=bytscl(tsel.halotimestep,top=254)
  loadct,0,/silent
  jeje_symbols,1,1.5
  plots,tsel.halotimestep,psel.mtot,psym=8,color=80
  
  return
  
end






