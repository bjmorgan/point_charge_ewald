!  errors.f90
!  28th November 2005
!  S. Reed 
!
! A simple subroutine to deal with errors in opening files.
! Gives an error message and stops the program

subroutine fileerror( ierr, filename )

  ! Get the error number and filename from the calling subroutine

  integer, intent( in )      :: ierr
  character(*), intent( in ) :: filename

  ! Give a specific message if we know what the error number means
  ! if not, give a message including the error number and stop.
 
  if (ierr == 10) then ! File already exists but was required to be new
     write(6,*) "Fatal error: The file ",trim(filename)," file already exists."
  else if (ierr == 29) then ! File doesn't exist but should
     write(6,*) "Fatal error: The file ",trim(filename)," does not exist so I can't open it!" 
  else ! Another error whose meaning we don't know 
     write(6,*) "Fatal error, ", ierr, ": Problem opening the file ", trim(filename)
  end if
  
  stop

end subroutine fileerror

subroutine filereaderror(ierr,errlocation,message)

  implicit none
  integer, intent(in) :: ierr
  character(*),intent(in) :: errlocation
  character(*), intent(in) :: message
  if (ierr == -1) then
     write(6,*) "Fatal error: Reached the end of file whilst reading from", trim(errlocation), "."
     stop
  else     
     write(6,*) "Fatal error number ",ierr," in ",trim(errlocation)
  end if

  write(6,*) "  ", trim(message)
  stop
end subroutine filereaderror
