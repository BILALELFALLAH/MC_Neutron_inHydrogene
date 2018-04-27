program trasport_scatering_abs_Hydrogene
implicit none
real::diff,diif
integer::i
real::Eint,sigmat,sigmad,sigmaa,phi,teta,r,L,x0,y0,z0,gsigmat,cm
   
x0=0;y0=0;z0=0 ; diif=0.
      do i=1,1000
          Eint=2.E6
          diff=0.
            do
             call g(Eint,sigmat,sigmad,sigmaa)
                 gsigmat=((6.023E23)*sigmat)/9  
                 CALL RANDOM_NUMBER(R)
                 teta=acos(2*R-1)
                 phi=2*4*atan(1.)*r
                 L=-log(R)/gsigmat
                 x0=x0+L*sin(teta)*cos(phi)            
                 y0=y0+L*sin(teta)*sin(phi)
                 z0=z0+L*cos(teta)
                 if (r<=(sigmaa/sigmat)) then
                    exit
                else    
                   cm=2*R-1      !cm=cos(tetacm)
                   Eint=(Eint*(1.+cm))/2
                    diff=diff+1
                    if(Eint<1.) then
                    exit
                    end if
                endif
          enddo
         diif=diif+diff
             print*,Eint,diff
    enddo
print*,diif/i
end program trasport_scatering_abs_Hydrogene
!---------------------------interpolation de section ef en fonction de energie---------------------------
subroutine g(Eint,sigmat,sigmad,sigmaa)
implicit none
real,dimension(307)::sigmatotal,Tn,sigmadiff,sigmaabso
real::sigmat,Eint,sigmad,sigmaa
integer::i
open(unit=2,file='SECTION_HY.txt')
do i=1,307
     read(2,*)Tn(i),sigmatotal(i),sigmadiff(i),sigmaabso(i)
enddo
close(2)
do i=1,307
    if(Eint>=Tn(i) .and. Eint<=Tn(i+1)) then
        sigmat=((sigmatotal(i)-sigmatotal(i+1))/(Tn(i)-Tn(i+1)))*(Eint-Tn(i))+sigmatotal(i)
        sigmad=((sigmadiff(i)-sigmadiff(i+1))/(Tn(i)-Tn(i+1)))*(Eint-Tn(i))+sigmadiff(i)
        sigmaa=((sigmaabso(i)-sigmaabso(i+1))/(Tn(i)-Tn(i+1)))*(Eint-Tn(i))+sigmaabso(i)
    endif
enddo
end subroutine 
	
	



