interface
   subroutine set_rule(mode, qmax, quad_idx, left, right, z, weight, zeta, lambda)
     implicit none
     integer, intent(in) :: mode, qmax, quad_idx
     COMPLEX_TYPE, intent(in) :: left
     REAL_TYPE, intent(in) :: right
     COMPLEX_TYPE, intent(out) :: z, weight, zeta
     COMPLEX_TYPE, intent(inout) :: lambda
   end subroutine set_rule
end interface
