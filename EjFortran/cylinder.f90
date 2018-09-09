program cylinder

! Calculate the surface area of a cylinder.
!
! Declare variables and constants.
! constants=pi
! variables=radius squared and height

  implicit none    ! Require all variables to be explicitly declared

  integer, parameter :: sp = selected_real_kind(6,37) 
  integer, parameter :: dp = selected_real_kind(15,307) 
  
  integer :: ierr, nSamples, idSample, i
  character(1) :: yn
  real(kind = dp) :: radius, height, area, promedio
  real(kind = dp), parameter :: pi = 4.0D0*atan(1.0D0)

  open(unit = 25, file = "datosF.ipn")
  read (25,*,iostat=ierr) nSamples

  interactive_loop: do i=1, nSamples

!   Prompt the user for radius and height
!   and read them.
    !write (*,*) 'Ingrese radio de la base y altura del cilindro: '

    read (25,*,iostat=ierr) idSample, radius, height


!   If radius and height could not be read from input,
!   then cycle through the loop.

    if (ierr /= 0) then
      write(*,*) 'Error, invalid input.'
      exit interactive_loop !Esto termina el loop
      !stop !esto detiene todo el programa
    end if

!   Compute area.  The ** means "raise to a power."

    area = 2.0D0*pi*(radius**2 + radius*height)

!   Write the input variables (radius, height)
!   and output (area) to the screen.

    write (*,'(1x,a12,i4,a9,f7.2,a9,f7.2,a7,f7.2)') &
      'Sample Id = ',idSample, ', radius = ', radius,', height = ',height,', area = ',area
    promedio = promedio + area
   

  end do interactive_loop

  promedio = promedio/dble(nSamples)
  write(*,"(a42,2x,f6.2)") "The average area for all the sample is = ",promedio
  close(25)
  stop
end program cylinder