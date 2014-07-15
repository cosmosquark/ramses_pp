PRO compare_branch


  file = 'progenitor.dat'
  read_progenitor,file,px,py,pz,zz,ns
  file = '/data/leo/CosmoZooms/Gorda-A2/progenitor.dat'
  read_progenitor,file,px2,py2,pz2,zz2,ns2

  ps = 0
  initcolor,ps=ps,eps=eps,black,blue,red,green,magenta,gc,gf,yellow,cyan,orange,purple
  window,retain=2,xsize=1000,ysize=1000
  !P.multi=[0,2,2]
  !P.charsize=2
  plot,zz,px,psym=4,xrange=[20.,0.1],/xlog,yrange=[0.4,0.6]
  oplot,zz2,px2,color=orange

  oplot,zz,py,psym=4
  oplot,zz2,py2,color=orange

  oplot,zz,pz,psym=4
  oplot,zz2,pz2,color=orange


  i1 = 6
  i2 = 184


  plot,ns[i1:i2],px[i1:i2]-px2,yrange=[-1.e-4,1.e-4]
  oplot,ns[i1:i2],py[i1:i2]-py2,color=green
  oplot,ns[i1:i2],pz[i1:i2]-pz2,color=cyan

  plot,px,py,psym=4,xrange=[0.4,0.6],yrange=[0.4,0.6]
  oplot,px2,py2,color=orange

  plot,px,pz,psym=4,xrange=[0.4,0.6],yrange=[0.4,0.6]
  oplot,px2,pz2,color=orange

  return
END

PRO read_progenitor,file,px,py,pz,zz,ns

  openr,uu,file,/get_lun
  nprog = 0L
  i1 = (i2 = (i3 = (i4 = 0L)))
  f1 = (f2 = (f3 = (f4 = (f5 = 0.0))))
  readf,uu,nprog
  px = (py = (pz = (zz = fltarr(nprog))))
  ns = intarr(nprog)
  for i=0,nprog-1 do begin
     ;;readf,uu,format='(3(2x,i5),2x,f8.5,2x,f6.3,2x,i8,3(2x,f9.6))',i1,i2,i3,f1,f2,i4,f3,f4,f5
     readf,uu,i1,i2,i3,f1,f2,i4,f3,f4,f5
     px(i) = f3
     py(i) = f4
     pz(i) = f5
     zz(i) = f2
     ns(i) = i2
  endfor
  free_lun,uu

  print,nprog
  print,i1,i2,i3,f1,f2,i4,f3,f4,f5

  return
END

