function read_halo,unit
  
  my_number=(my_timestep=(level=(hosthalo=(hostsub=(nbsub=(nextsub=0L))))))
  m=(px=(py=(pz=(vx=(vy=(vz=(Lx=(Ly=(Lz=(rr=(a=(b=(c=(ek=(ep=(et=(spin=(rvir=(mvir=(tvir=(cvel=(rho0=(r_c=0.0)))))))))))))))))))))))
  r = 0L
  readu,unit,r,my_number,r
  readu,unit,r,my_timestep,r
  readu,unit,r,level,hosthalo,hostsub,nbsub,nextsub,re
  readu,unit,r,m,r
  readu,unit,r,px,py,pz,r
  readu,unit,r,vx,vy,vz,r
  readu,unit,r,lx,ly,lz,r
  readu,unit,r,rr,a,b,c,r
  readu,unit,r,ek,ep,et,r
  readu,unit,r,spin,r
  readu,unit,r,rvir,mvir,tvir,cvel,r
  readu,unit,r,rho0,r_c,r

  h = {my_number:my_number,my_timestep:my_timestep,level:level,hosthalo:hosthalo,hostsub:hostsub,nbsub:nbsub,nextsub:nextsub, $
       m:m,px:px,py:py,pz:pz,vx:vx,vy:vy,vz:vz,Lx:lx,Ly:ly,Lz:lz, $
       r:rr,a:a,b:b,c:c,ek:ek,ep:ep,et:et,spin:spin,rvir:rvir,mvir:mvir,tvir:tvir,cvel:cvel,rho0:rho0,r_c:r_c}
  
  return,h

end 

pro read_treebrick,file,treeb,minPartID,swap=swap,stars=stars

openr,11,file
r = (npart = (nhalo = (nsub = 0L)))
x = 0.0
readu,11,r,npart,r 
readu,11,r,x,r
readu,11,r,x,r
readu,11,r,x,r
readu,11,r,x,r
readu,11,r,nhalo,nsub,r 
if keyword_set(swap) then begin
   npart = swap_endian(npart)
   nhalo = swap_endian(nhalo)
   nsub  = swap_endian(nsub)
endif
print,npart,nhalo,nsub

halo = {my_number:0L,my_timestep:0L,level:0L,hosthalo:0L,hostsub:0L,nbsub:0L,nextsub:0L, $
       m:0.0,px:0.0,py:0.0,pz:0.0,vx:0.0,vy:0.0,vz:0.0,Lx:0.0,Ly:0.0,Lz:0.0, $
       r:0.0,a:0.0,b:0.0,c:0.0,ek:0.0,ep:0.0,et:0.0,spin:0.0,rvir:0.0,mvir:0.0,tvir:0.0, $
       cvel:0.0,rho0:0.0,r_c:0.0}
treeb = replicate(halo,nhalo+nsub)
minPartID = lon64arr(nhalo+nsub)

np = 0L
ih = 0L
while not eof(11) do begin 
   readu,11,r,np,r 
   if keyword_set(swap) then np = swap_endian(np)
   listp = lonarr(np) 
   readu,11,r,listp,r
   if (keyword_set(swap)) then listp=swap_endian(listp)
   minPartID(ih) = min(listp)
   if keyword_set(stars) then begin 
      ;; read stars information
      a = fltarr(np)
      readu,11,r,a,r  ;; mass
      print,'mass: ',min(a),max(a)
      readu,11,r,a,r  ;; age (in fact conformal time...)
      print,'age: ',min(a),max(a)
      readu,11,r,a,r  ;; metallicity
      print,'met: ',min(a),max(a)
   endif
   treeb(ih) = read_halo(11)
   if keyword_set(swap) then treeb(ih) = swap_endian(treeb(ih))
   ih++
endwhile

close,11

end



