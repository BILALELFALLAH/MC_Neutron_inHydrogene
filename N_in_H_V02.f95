! by bilal el fallah 
program N_in_H
implicit none 
real ::siigmma_scatering,siigmma_absorption,siigmma_total
real ::nbr_scatering,nbr_scatering_moy,E_int
real :: COS_THETA_CM ,R 
integer ::i
print*,'   E_out         Nbr_scatering' 
nbr_scatering_moy=0.
do i=1,500
         E_int=2.E6 
         nbr_scatering=0.
         do
 CALL g(E_int,siigmma_absorption,siigmma_scatering,siigmma_total)
               CALL RANDOM_NUMBER(R)
       if(R <= siigmma_absorption/siigmma_total) then 
       exit
       else        
               COS_THETA_CM=-1+2*R
               E_int=E_int*(1+COS_THETA_CM)/2
               nbr_scatering=nbr_scatering+1
               if (E_int <= 1.) then 
               exit
               end if
               endif
         enddo
    print*,E_int,nbr_scatering  
   nbr_scatering_moy = nbr_scatering_moy + nbr_scatering
end do
          print*,'the mean number of scatering is :',nbr_scatering_moy/i                    
end program N_in_H
subroutine  g(E_int,siigmma_absorption,siigmma_scatering,siigmma_total)
implicit none 
real,dimension(307):: sigmma_absorption,sigmma_scatering,sigmma_total,Tn
real  :: siigmma_absorption,siigmma_scatering,siigmma_total,E_int
integer :: i
open(unit=4,file='SECTION_HY.txt')
do i=1,307
read(4,*)Tn(i),sigmma_total(i),sigmma_scatering(i),sigmma_absorption(i)
end do 
close(4)
     do i=1,307
            if (E_int>=Tn(i) .and. E_int <= Tn(i+1))  then 
siigmma_total=((sigmma_total(i)-sigmma_total(i+1))/(Tn(i)-Tn(i+1)))*(E_int-Tn(i))+sigmma_total(i)
siigmma_scatering=((sigmma_scatering(i)-sigmma_scatering(i+1))/(Tn(i)-Tn(i+1)))*(E_int-Tn(i))+sigmma_scatering(i)
siigmma_absorption=((sigmma_absorption(i)-sigmma_absorption(i+1))/(Tn(i)-Tn(i+1)))*(E_int-Tn(i))+sigmma_absorption(i)              
           endif
           end do
end subroutine           


          
