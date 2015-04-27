!
!  Binary tree cache of nuclear evolution coefficients a-la 
!  Gruner and Brumer (see franck_condon.f90 for the appropriate reference).
!
!  The integrals are stored in a binary-tree structure, with the
!  path from the root of the tree determining the state index of
!  a given node. The tree could only grow - so we better have 
!  enough physical memory to accommodate the entire structure.
!
!  Note that cache_gr_initialize and cache_gr_destroy MUST be called from
!  a serial serion. cache_gr_locate may be called from a parallel region,
!  but all updates will be serialized.
!
!  WARNING: This module is only thread-safe on cache-coherent architectures.
!
  module cache_gr
    use accuracy
    use timer
    !$ use OMP_LIB
    implicit none
    private
    public grCache, grCacheWrapper
    public cache_gr_invalid_r, cache_gr_invalid_c
    public cache_gr_initialize, cache_gr_destroy, cache_gr_locate
    !
    real(rk), parameter    :: cache_gr_invalid_r = huge(1._rk)
    complex(rk), parameter :: cache_gr_invalid_c = (cache_gr_invalid_r,cache_gr_invalid_r)
    type grCache
      complex(rk), allocatable  :: values(:)       ! Taylor series coefficients, starting at
                                                   ! zero-th order. Uninitialized values are set
                                                   ! to cache_gr_invalid.
      type(grCache), pointer    :: left  => NULL() ! Increment current excitation level
      type(grCache), pointer    :: right => NULL() ! Increment mode number
      !$ integer(OMP_LOCK_KIND) :: update_lock     ! Must be locked before allocation/deallocation 
                                                   ! takes place
    end type ! grCache
    !
    type grCacheWrapper
      type(grCache), pointer         :: root => NULL()   ! Tree of cached FC integrals
      integer(ik)                    :: order            ! Maximum order of Taylor coefficients
      integer(sik), allocatable      :: n_max(:)         ! Max. allowed excitation order in the index
    end type ! grCacheWrapper
    !
    contains
    !
    subroutine cache_gr_initialize(cache,order,n_max)
      type(grCacheWrapper), intent(inout) :: cache     ! Cache tree descriptor
      integer(ik), intent(in)             :: order     ! Max. order of Taylor coefficients
      integer(sik), intent(in)            :: n_max(:)  ! Max. allowed excitation orders in the index
      !
      integer(ik) :: alloc1, alloc2
      !
      if (associated(cache%root)) then
        stop 'cache_gr%cache_gr_initialize - double initialization'
      end if
      !
      allocate (cache%root,cache%n_max(size(n_max)),stat=alloc1) 
      allocate (cache%root%values(0:order),stat=alloc2)
      if (alloc1/=0 .or. alloc2/=0) then
        write (out,"('Error ',i8,'/',i8,' allocating coefficient cache root')") alloc1, alloc2
        stop 'cache_gr%cache_gr_initialize - allocation failed'
      end if
      cache%n_max = n_max
      cache%order = order
      cache%root%left   => NULL()
      cache%root%right  => NULL()
      cache%root%values = cache_gr_invalid_c
      !$ call omp_init_lock(cache%root%update_lock)
    end subroutine cache_gr_initialize
    !
    subroutine cache_gr_destroy(cache,tree_size)
      type (grCacheWrapper), intent(inout) :: cache      ! Cache descriptor
      real(ark), intent(out)               :: tree_size  ! Number of nodes in the tree
      !
      integer(ik) :: alloc
      !
      call gr_destroy_tree(cache%root,tree_size)
      deallocate (cache%n_max,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' deallocating max. index arrays')") alloc
        stop 'cache_gr%cache_gr_destroy - deallocation failed'
      end if
    end subroutine cache_gr_destroy
    !
    recursive subroutine gr_destroy_tree(root,tree_size)
      type(grCache), pointer :: root       ! Tree to be destroyed
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
      call gr_destroy_tree(root%left,left_size)
      call gr_destroy_tree(root%right,right_size)
      if (allocated(root%values)) then
        deallocate (root%values,stat=alloc)
      end if
      !$ call omp_destroy_lock(root%update_lock)
      deallocate (root,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' deallocating the tree')") alloc
        stop 'cache_gr%gr_destroy_tree - error deallocating the tree'
      end if
      root => NULL()
      tree_size = 1 + left_size + right_size
    end subroutine gr_destroy_tree
    !
    function cache_gr_locate(cache,n) result(cgr)
      type(grCacheWrapper), intent(inout) :: cache  ! Cache descriptor
      integer(sik), intent(in)            :: n(:)   ! Indices for the amplitude
      type(grCache), pointer              :: cgr    ! Cache entry. If any of the excitation indices
                                                    ! exceed the allowed range, cfc will be set to NULL()
      !
      integer(sik) :: nc(size(cache%n_max)) ! Packed index
      integer(ik)  :: n_top
      integer(ik)  :: n_end                 ! The last non-zero element of nc
      !
      !  Is the requested amplitude cacheable? If not, return failure
      !
      if (any(n<0) .or. any(n>cache%n_max) ) then
        cgr => NULL()
        return
      end if
      !
      !  Pack significant state indices together
      !
      n_top = count(cache%n_max>0)
      nc(:n_top) = pack(n,mask=(cache%n_max>0))
      !
      find_end: do n_end=n_top,1,-1
        if (nc(n_end)>0) exit find_end
      end do find_end
      ! write(out,"('full: ',40i3)") n
      ! write(out,"('pack: ',40i3)") nc(:n_end)
      !
      !  Look up the amplitude, allocating missing nodes along the way
      !
      cgr => gr_locate(cache%root,nc(:n_end))
    end function cache_gr_locate
    !
    !  We need to be a little careful with parallel updates - there is 
    !  a possibility of a race condition between association test and 
    !  acquisition of the lock - so all association tests have to be 
    !  repeated after acquiring the update lock. Furthermore, we cannot
    !  connect the new node to the cache tree until we have finished
    !  it's initialization in toto, since we do not prohibit lookups
    !  of the node being updated.
    !
    recursive subroutine gr_grow_node_right(node)
      type(grCache), pointer :: node  ! Node to update
      !
      integer(ik)            :: lb, ub, alloc1, alloc2
      type(grCache), pointer :: tmp
      !
      !$ call omp_set_lock(node%update_lock)
      !$ if (.not.associated(node%right)) then
        lb = lbound(node%values,dim=1)
        ub = ubound(node%values,dim=1)
        allocate (tmp,stat=alloc1)
        allocate (tmp%values(lb:ub),stat=alloc2)
        if (alloc1/=0 .or. alloc2/=0) then
          write (out,"('Error growing cache tree to the right: ',i8,1x,i8)") alloc1, alloc2
          stop 'cache_gr%gr_grow_node_right - allocation failed'
        end if
        tmp%left   => NULL()
        tmp%right  => NULL()
        tmp%values(:) = cache_gr_invalid_c
        !$ call omp_init_lock(tmp%update_lock)
        !
        !  The node is fully initialized; it is safe to connect it even before we
        !  release the lock. The assignment _must_ be atomic, or bad things 
        !  will happen. Unfortunately, omp atomic does not work for pointer
        !  assignments - so we'll just have to pray.
        !
        node%right => tmp
      !$ end if
      !$ call omp_unset_lock(node%update_lock)
    end subroutine gr_grow_node_right
    !
    recursive subroutine gr_grow_node_left(node)
      type(grCache), pointer :: node  ! Node to update
      !
      integer(ik)            :: lb, ub, alloc1, alloc2
      type(grCache), pointer :: tmp
      !
      !$ call omp_set_lock(node%update_lock)
      !$ if (.not.associated(node%left)) then
        lb = lbound(node%values,dim=1)
        ub = ubound(node%values,dim=1)
        allocate (tmp,stat=alloc1)
        allocate (tmp%values(lb:ub),stat=alloc2)
        if (alloc1/=0 .or. alloc2/=0) then
          write (out,"('Error growing cache tree to the left: ',i8,1x,i8)") alloc1, alloc2
          stop 'cache_gr%gr_grow_node_left - allocation failed'
        end if
        tmp%left   => NULL()
        tmp%right  => NULL()
        tmp%values(:) = cache_gr_invalid_c
        !$ call omp_init_lock(tmp%update_lock)
        node%left => tmp
      !$ end if
      !$ call omp_unset_lock(node%update_lock)
    end subroutine gr_grow_node_left
    !
    recursive function gr_locate(node,nc) result(cgr)
      type(grCache), pointer      :: node  ! Current node
      integer(sik), intent(inout) :: nc(:) ! Current state vector
      type(grCache), pointer      :: cgr
      !
      if (size(nc)<=0) then
        !
        !  This is the mode we've been looking for
        !
        cgr => node
      else if (nc(1)<=0) then
        !
        !  We have enough quanta in the current mode. Are there more modes?
        !
        if (size(nc)<=1) then
          !
          !  No - this is the mode we've been looking for
          !
          cgr => node
        else
          !
          !  We have more modes - recurse to the right.
          !
          if (.not.associated(node%right)) then
            call gr_grow_node_right(node)
          end if
          cgr => gr_locate(node%right,nc(2:))
        end if
      else ! if (nc(1)<=0)
        !
        !  Need more quanta in the current mode - recurse to the left
        !
        nc(1) = nc(1) - 1
        if (.not.associated(node%left)) then
          call gr_grow_node_left(node)
        end if
        cgr => gr_locate(node%left,nc)
        nc(1) = nc(1) + 1
      end if
      !
    end function gr_locate
  end module cache_gr
