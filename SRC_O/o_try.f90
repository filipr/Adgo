module try_mod
   use time_mod

   contains

   subroutine try(this)
         class (Time_t), intent( inout ) :: this

         print*, 'try'
   end subroutine try

end module try_mod
