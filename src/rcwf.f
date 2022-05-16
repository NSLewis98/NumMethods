subroutine RCWF(rho,eta,l,fc,gc,fcp,gcp)
  !
  !.. Description
  !..
  ! Coulomb wavefunctions calculated at r = rho by the
  ! continued-fraction method of Steed. See:
  ! A. R. Barnett et al, Comput. Phys. Commun.  8, 377(1974) 
  ! A. R. Barnett      , Comput. Phys. Commun. 11, 141(1976) 
  !
  implicit none
  !
  !.. Input variables
  !..
  real(kind(1d0)), intent(in) :: rho,eta
  real(kind(1d0)), intent(out):: fc,gc,fcp,gcp
  integer        , intent(in) :: l
  !
  !.. Parameters
  !..
  real(kind(1d0)), parameter :: ACCURACY =   1.d-15
  real(kind(1d0)), parameter :: STEP     = 999.d0
  !
  !.. Local variables
  !..
  real(kind(1d0)) :: xll1,eta2,turn,ktrp,pmx,fp,rh
  real(kind(1d0)) :: rho2,pl,tf,dk,del,d,h,w,g,gp,f
  real(kind(1d0)) :: r3,p,q,ar,ai,br,bi,wi,dr,di
  real(kind(1d0)) :: dp,dq,t,s,h2,i2,etah,h2ll,rh2
  real(kind(1d0)) :: r,ktr,etar,tfp
  real(kind(1d0)) :: k,k1,k2,k3,k4,m1,m2,m3,m4
  integer         :: iflag

  r    = rho                                                        
  ktr  = 1                                                          
  xll1 = dble(l*(l+1))                                          
  eta2 = eta*eta                                                    
  turn = eta + sqrt(eta2 + xll1)                                    
  if(r.lt.turn.and.abs(eta).ge.1.0e-6) ktr = -1                     
  ktrp = ktr                                                        

  iflag=0
  outer : do
     do
        if(iflag==1)then
           r    = turn                                                       
           tf   = f                                                          
           tfp  = fp                                                         
           ktrp = 1     
        endif
        iflag=1
        do
           etar = eta*r                                                      
           rho2 =   r*r                                                      
           pl   = dble(l+1)                                            
           pmx  = pl + 0.5                                                   
           !continued fraction for fp(l)/f(l)
           !xl is f  xlprime is fp **
           fp  = eta/pl + pl/r                                               
           dk  = etar*2.0                                                    
           del = 0.0                                                         
           d   = 0.0                                                         
           f   = 1.0                                                         
           k   = (pl*pl - pl + etar)*(2.0*pl - 1.0)                          
           if(abs(pl*pl+pl+etar)>ACCURACY)exit
           r   = r + 1.0e-6                                                  
        enddo
        do
           h   = (pl*pl + eta2)*(1.0 - pl*pl)*rho2                           
           k   = k + dk + pl*pl*6.0                                          
           d   =  1.0/(d*h + k)                                              
           del =  del*(d*k - 1.0)                                            
           if(pl.lt.pmx) del = -r*(pl*pl + eta2)*(pl + 1.0)*d/pl             
           pl  = pl + 1.0                                                    
           fp  = fp + del                                                    
           if(d.lt.0.0) f = -f                                               
           if(pl.gt.20000.)then
              w  = 0.0
              g  = 0.0
              gp = 0.0
              exit outer
           endif

           if(abs(del)<abs(fp)*ACCURACY)exit
        enddo
        fp  = f*fp                                                        
        if(abs(ktrp+1.d0)>ACCURACY)exit
     enddo

     !repeat for r = turn if rho lt turn                                
     !now obtain p + i.q for l from continued fraction (32)          
     !real arithmetic to facilitate conversion to ibm using real*8      
     p  = 0.0                                                          
     q  = r - eta                                                      
     pl = 0.0                                                          
     ar = -(eta2 + xll1)                                               
     ai =   eta                                                        
     br = 2.0*q                                                        
     bi = 2.0                                                          
     wi = 2.0*eta                                                      
     dr =   br/(br*br + bi*bi)                                         
     di =  -bi/(br*br + bi*bi)                                         
     dp = -(ar*di + ai*dr)                                             
     dq =  (ar*dr - ai*di)                                             
     do
        p  =  p + dp                                                      
        q  =  q + dq                                                      
        pl = pl + 2.0                                                     
        ar = ar + pl                                                      
        ai = ai + wi                                                      
        bi = bi + 2.0                                                     
        d  = ar*dr - ai*di + br                                           
        di = ai*dr + ar*di + bi                                           
        t  = 1.0/(d*d + di*di)                                            
        dr =  t*d                                                         
        di = -t*di                                                        
        h  = br*dr - bi*di - 1.0                                          
        k  = bi*dr + br*di                                                
        t  = dp*h  - dq*k                                                 
        dq = dp*k  + dq*h                                                 
        dp = t                                                            
        if(rho==0.d0.and.pl>46000.)then
           w  = 0.0                                                          
           g  = 0.0                                                          
           gp = 0.0                                                          
           exit outer
        endif
        if(abs(dp)+abs(dq)<(abs(p)+abs(q))*ACCURACY)exit
     enddo
     p  = p/r                                                          
     q  = q/r                                                          
     !solve for fp,g,gp and normalise f  at l=l                      
     g  = (fp - p*f)/q                                                 
     gp = p*g - q*f                                                    
     w  = 1.0/sqrt(fp*g - f*gp)                                        
     g  = w*g                                                          
     gp = w*gp                                                         
     if(abs(ktr-1.d0)<ACCURACY)exit outer
     f  = tf                                                           
     fp = tfp                                                          
     !runge-kutta integration of g(l) and gp(l) inwards from turn 
     !            see fox and mayers 1968 pg 202                        
     r3 = 1.0/3.0d0                                                    
     h  = (rho - turn)/(STEP + 1.0)                                    
     h2 = 0.5d0*h                                                      
     i2 = int(STEP + 0.001)                                           
     etah = eta*h                                                      
     h2ll = h2*xll1                                                    
     s  = (etah + h2ll/r  )/r   - h2                                   
     do
        rh2= r + h2                                                       
        t  = (etah + h2ll/rh2)/rh2 - h2                                   
        k1 = h2*gp                                                        
        m1 =  s*g                                                         
        k2 = h2*(gp + m1)                                                 
        m2 =  t*(g  + k1)                                                 
        k3 =  h*(gp + m2)                                                 
        m3 =  t*(g  + k2)                                                 
        m3 =     m3 + m3                                                  
        k4 = h2*(gp + m3)                                                 
        rh = r + h                                                        
        s  = (etah + h2ll/rh )/rh  - h2                                   
        m4 =  s*(g + k3)                                                  
        g  = g  + (k1 + k2 + k2 + k3 + k4)*r3                             
        gp = gp + (m1 + m2 + m2 + m3 + m4)*r3                             
        r  = rh                                                           
        i2 = i2 - 1                                                       
        if(abs(gp).gt.1.d300)then
           w  = 0.0                                                          
           g  = 0.0                                                          
           gp = 0.0                                                          
           exit outer
        endif
        if(i2<0)exit
     enddo
     w  = 1.0/(fp*g - f*gp)                                            
     exit outer
  enddo outer

  !gc and gcp, stored values are r,s
  !renormalise fc,fcp                            
  gc =g                                       
  gcp=gp                                      
  fcp=fp*w                                    
  fc= f*w

  return                                                            

end subroutine RCWF
