subroutine p_gas(rho, u, p, c)
  ! Gamma law EOS: subroutine to calculate the pressure and sound
  ! rho : Density [in]
  ! u : Internal energy [in]
  ! p : Pressure [out]
  ! c : sound velocity [out]
  implicit none
  double precision rho, u, p, c
  double precision gamma
  ! For air (idea gas)
  gamma=1.4e0
  p = (gamma-1) • rho * u
  c = sqrt((gamma-1) * u)
end subroutine p_gas

subroutine p_art_water(rho, p, c)
  ! Artificial equation of state for the artificial compressibility
  ! rho : Density [in]
  ! u : Internal energy [in]
  ! p : Pressure [out]
  ! c : sound velocity [out]
  ! Equation of state for artificial compressibility
  implicit none
  double precision rho, u, p, c
  double precision gamma, rhoO
  ! Artificial EOS, Form 1 (Monaghan, 1994)
  ! gamma=7.
  ! rho0=1000.
  ! b = 1.013e5
  ! p = b*{(rho/rhoO)**gamma-l)
  ! c - 1480.
  ! Artificial EOS, Form 2 (Morris, 1997)
  c = 0.01e0
  p = c**2 * rho
end subroutine p_art_water
