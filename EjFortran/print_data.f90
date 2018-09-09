program data
 implicit none

 integer, parameter :: sp = selected_real_kind(6,37)
 integer, parameter :: dp = selected_real_kind(15,307)
 integer :: i, nSamples

 real :: rprom, hprom, A1, A2
 real(kind=dp) :: r,h
 real(kind=8) :: n

 write(*,*) "Ingrese el radio promedio ,altura promedio y numero de datos"
 read(*,*) rprom,hprom,nSamples

 A1 = 1.0; A2 = 1.0; r = 5
 
 open(10, file = "datosF.ipn")
 escritura: do i = 1, nSamples
  
  r = rprom + A1*sin(dble(i)/10); h = rprom + A2*cos(dble(i)/10)
  call random_number(n)
  write(10,*)int(100*n),r,h

 end do escritura
 write(*,*) "Done!"
end program data