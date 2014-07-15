pro plot_accretion_rate,timestep,datadir

  filename = strtrim(datadir,2) + '/halos_results.' + string(timestep,format='(i3.3)')
  r=(ts=(nh=(ns=(nsteps=0L))))
  aexp=0.0
  openr,11,filename
  readu,11,r,ts,nh,ns,aexp,nsteps,r
  
  hst = {level:0L,hosthalo:0L,hostsub:0L,nbsub:0L,nextsub:0L,x:0.0,y:0.0,z:0.0, $
         vx:0.0,vy:0.0,vz:0.0,m:0.0,r:0.0,spin:0.0,sha:0.0,shb:0.0,shc:0.0,     $
         et:0.0,ek:0.0,ep:0.0,lx:0.0,ly:0.0,lz:0.0,rvir:0.0,mvir:0.0,tvir:0.0,  $
         cvel:0.0,rho:0.0,rc:0.0,macc:0.0d0}
  h   = replicate(hst,nh+ns)
  readu,11,r,h,r
  close,11

  stop
end
