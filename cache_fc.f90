!
!  Binary tree cache of Franck-Condon integrals a-la Gruner and Brumer
!  (see franck_condon.f90 for the appropriate reference).
!
!  The integrals are stored in a binary-tree structure, with the
!  path from the root of the tree determining the state index of
!  a given node. The tree could only grow - so we better have 
!  enough physical memory to accommodate the entire structure.
!
!  WARNING: This module is not thread-safe!
!
  module cache_fc
    use accuracy
    use timer
    implicit none
    private
    public fcCache, fcCacheWrapper
    public cache_fc_invalid
    public cache_fc_initialize, cache_fc_destroy, cache_fc_locate
    !
    real(rk), parameter :: cache_fc_invalid = 100._rk
    type fcCache
      real(rk)               :: value = cache_fc_invalid ! Value of the FC integral. 
      real(rk)               :: weight = -1._rk          ! Significance level of the cached value
      type(fcCache), pointer :: left  => NULL()          ! Increment current excitation level
      type(fcCache), pointer :: right => NULL()          ! Increment mode number
    end type ! fcCache
    !
    type fcCacheWrapper
      type(fcCache), pointer    :: root => NULL()   ! Tree of cached FC integrals
      integer(ik)               :: state_size       ! Total number of state indices
      integer(ik)               :: n2_size, n1_size ! ... and their n2/n1 components
      integer(sik), allocatable :: n2_max(:)        ! Max. allowed order in lhs
      integer(sik), allocatable :: n1_max(:)        ! Max. allowed order in rhs
    end type ! fcCacheWrapper
    !
    contains
    !
    subroutine cache_fc_initialize(cache,n2_max,n1_max)
      type(fcCacheWrapper), intent(inout) :: cache     ! Cache tree descriptor
      integer(sik), intent(in)            :: n2_max(:) ! Max. allowed indices in lhs
      integer(sik), intent(in)            :: n1_max(:) ! Max. allowed indices in rhs
      !
      integer(ik) :: alloc
      !
      if (associated(cache%root)) then
        stop 'cache_fc%cache_fc_initialize - double initialization'
      end if
      !
      allocate (cache%root,cache%n2_max(size(n2_max)),cache%n1_max(size(n1_max)),stat=alloc) 
      if (alloc/=0) then
        write (out,"('Error ',i8,' allocating FC cache root')") alloc
        stop 'cache_fc%cache_fc_initialize - allocation failed'
      end if
      cache%n2_max = n2_max
      cache%n1_max = n1_max
      cache%n2_size = count(n2_max>0)
      cache%n1_size = count(n1_max>0)
      cache%state_size = cache%n2_size + cache%n1_size
      cache%root%value = cache_fc_invalid
      cache%root%weight= -1._rk
      cache%root%left => NULL()
      cache%root%right=> NULL()
    end subroutine cache_fc_initialize
    !
    subroutine cache_fc_destroy(cache,tree_size)
      type (fcCacheWrapper), intent(inout) :: cache      ! Cache descriptor
      real(ark), intent(out)               :: tree_size  ! Number of nodes in the tree
      !
      integer(ik) :: alloc
      !
      call fc_destroy_tree(cache%root,tree_size)
      deallocate (cache%n2_max,cache%n1_max,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' deallocating max. index arrays')") alloc
        stop 'cache_fc%cache_fc_destroy - deallocation failed'
      end if
    end subroutine cache_fc_destroy
    !
    recursive subroutine fc_destroy_tree(root,tree_size)
      type(fcCache), pointer :: root       ! Tree to be destroyed
      real(ark), intent(out) :: tree_size  ! Number of nodes in the tree
      !
      real(ark)   :: left_size, right_size
      integer(ik) :: alloc
      !
      if (.not.associated(root)) then
        tree_size = 0
        return
      end if
      !
      call fc_destroy_tree(root%left,left_size)
      call fc_destroy_tree(root%right,right_size)
      deallocate (root,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' deallocating the tree')") alloc
        stop 'cache_fc%fc_destroy_tree - error deallocating the tree'
      end if
      root => NULL()
      tree_size = 1 + left_size + right_size
    end subroutine fc_destroy_tree
    !
    function cache_fc_locate(cache,n2,n1) result(cfc)
      type(fcCacheWrapper), intent(inout) :: cache  ! Cache descriptor
      integer(sik), intent(in)            :: n2(:)  ! Indices for the LHS
      integer(sik), intent(in)            :: n1(:)  ! and RHS
      type(fcCache), pointer              :: cfc    ! Cache entry. If any of the n2/n1 indices
                                                    ! exceed the allowed range, cfc will be set to NULL()
      !
      integer(sik) :: nc(cache%state_size)
      integer(ik)  :: n_end                 ! The last non-zero element of nc
      !
      !  Is the requested integral cacheable? If not, return failure
      !
      if (any(n2<0) .or. any(n2>cache%n2_max) .or.  any(n1<0) .or. any(n1>cache%n1_max)) then
        cfc => NULL()
        return
      end if
      !
      !  Pack significant state indices together
      !
      nc(1:cache%n2_size)  = pack(n2,mask=(cache%n2_max>0))
      nc(cache%n2_size+1:) = pack(n1,mask=(cache%n1_max>0))
      !
      find_end: do n_end=size(nc),1,-1
        if (nc(n_end)>0) exit find_end
      end do find_end
      !
      !  Look up the combined integral, allocating missing nodes along the way
      !
      cfc => fc_locate(cache%root,nc(:n_end))
    end function cache_fc_locate
    !
    recursive function fc_locate(node,nc) result(cfc)
      type(fcCache), pointer      :: node  ! Current node
      integer(sik), intent(inout) :: nc(:) ! Current state vector
      type(fcCache), pointer      :: cfc
      !
      if (size(nc)<=0) then
        !
        !  This is the mode we've been looking for
        !
        cfc => node
      else if (nc(1)<=0) then
        !
        !  We have enough quanta in the current mode. Are there more modes?
        !
        if (size(nc)<=1) then
          !
          !  No - this is the mode we've been looking for
          !
          cfc => node
        else
          !
          !  We have more modes - recurse to the right
          !
          if (.not.associated(node%right)) then
            allocate (node%right)
            node%right%value = cache_fc_invalid
            node%right%weight= -1._rk
            node%right%left  => NULL()
            node%right%right => NULL()
          end if
          cfc => fc_locate(node%right,nc(2:))
        end if
      else ! if (nc(1)<=0)
        !
        !  Need more quanta in the current mode - recurse to the left
        !
        nc(1) = nc(1) - 1
        if (.not.associated(node%left)) then
          allocate (node%left)
          node%left%value = cache_fc_invalid
          node%left%weight= -1._rk
          node%left%left  => NULL()
          node%left%right => NULL()
        end if
        cfc => fc_locate(node%left,nc)
        nc(1) = nc(1) + 1
      end if
    end function fc_locate
  end module cache_fc
