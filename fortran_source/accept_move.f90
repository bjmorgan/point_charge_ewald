logical function accept_move( new_energy, old_energy, kT )

implicit none

double precision, intent(in) :: new_energy, old_energy, kT
double precision :: boltzmann_factor, delta_E, r

delta_E = new_energy - old_energy
if ( delta_E .le. 0.0 ) then
    accept_move = .true.
else
    boltzmann_factor = exp( -delta_E / kT )
    if ( rand(0) .le. boltzmann_factor ) then
        accept_move = .true.
    else
        accept_move = .false.
    endif
endif

end function accept_move
