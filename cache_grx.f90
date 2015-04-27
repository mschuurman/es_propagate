!
!  This module is a drop-in extension of cache_gr.f90; it adds a tree walker.
!  It is safe to grow the tree during the walk; however, some of the new nodes
!  may be missed during the walk-through. It is safe to call walker_gr_next()
!  from a parallel region; however, this is achieved by serializing it internally.
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
  module cache_grx
    use accuracy
    use timer
    !$ use OMP_LIB
    implicit none
    private
    public grCache, grCacheWrapper, grCacheWalker
    public cache_gr_invalid_r, cache_gr_invalid_c
    public cache_gr_initialize, cache_gr_destroy, cache_gr_locate
    public walker_gr_initialize, walker_gr_destroy, walker_gr_next
    public walker_complete ! Debug/Example - enumerate all nodes in the cache
    !
    real(rk), parameter    :: cache_gr_invalid_r = huge(1._rk)
    complex(rk), parameter :: cache_gr_invalid_c = (cache_gr_invalid_r,cache_gr_invalid_r)
    type grCache
      complex(rk), allocatable  :: values(:,:)     ! expansion coefficients for current, previous
                                                   ! and updated time steps. Uninitialized values are set
                                                   ! to cache_gr_invalid.
      type(grCache), pointer    :: left  => NULL() ! Increment current excitation level
      type(grCache), pointer    :: right => NULL() ! Increment mode number
      !$ integer(OMP_LOCK_KIND) :: update_lock     ! Must be locked before allocation/deallocation 
                                                   ! takes place
    end type ! grCacheA
    !
    type pgrCache
      type(grCache), pointer    :: p               ! Fortran does not make it easy to declare arrays of pointers ...
    end type pgrCache
    !
    type grCacheWrapper
      type(grCache), pointer         :: root => NULL()   ! Tree of cached FC integrals
      integer(ik)                    :: order            ! Maximum order of Taylor coefficients
      integer(sik), allocatable      :: n_max(:)         ! Max. allowed excitation order in the index
    end type ! grCacheWrapper
    !
    type grCacheWalker
      type(grCacheWrapper), pointer  :: cache => NULL()  ! Tree structure we are walking through
      type(grCache), pointer         :: node             ! Node we inspected on the previous iteration
      type(pgrCache), allocatable    :: n_stack(:)       ! Stack of nodes
      integer(ik), allocatable       :: q_stack(:)       ! Stack of mode indices
      integer(ik)                    :: p_stack          ! Current position within the stack
      logical                        :: done             ! Walk is complete
      integer(sik), allocatable      :: quanta(:)        ! Number of quanta in the current node
      integer(ik)                    :: p_quanta         ! Mode currently being updated; modes with zero limit
                                                         ! on the number of quanta is are excluded
      !$ integer(OMP_LOCK_KIND)      :: traverse_lock    ! This lock must be acquired while we are walking the tree
    end type ! grCacheWalker
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
      allocate (cache%root%values(0:order,2),stat=alloc2)
      if (alloc1/=0 .or. alloc2/=0) then
        write (out,"('Error ',i8,'/',i8,' allocating coefficient cache root')") alloc1, alloc2
        stop 'cache_gr%cache_gr_initialize - allocation failed'
      end if
      cache%n_max = n_max
      cache%order = order
      cache%root%left   => NULL()
      cache%root%right  => NULL()
      cache%root%values = (0._rk,0._rk)
!      cache%root%values = cache_gr_invalid_c
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
        allocate (tmp%values(lb:ub,2),stat=alloc2)
        if (alloc1/=0 .or. alloc2/=0) then
          write (out,"('Error growing cache tree to the right: ',i8,1x,i8)") alloc1, alloc2
          stop 'cache_gr%gr_grow_node_right - allocation failed'
        end if
        tmp%left   => NULL()
        tmp%right  => NULL()
        tmp%values = (0._rk,0._rk)
!        tmp%values = cache_gr_invalid_c
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
        allocate (tmp%values(lb:ub,2),stat=alloc2)
        if (alloc1/=0 .or. alloc2/=0) then
          write (out,"('Error growing cache tree to the left: ',i8,1x,i8)") alloc1, alloc2
          stop 'cache_gr%gr_grow_node_left - allocation failed'
        end if
        tmp%left   => NULL()
        tmp%right  => NULL()
        tmp%values = (0._rk,0._rk)
!        tmp%values = cache_gr_invalid_c
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
    !
    subroutine walker_gr_initialize(cache,walker)
      type(grCacheWrapper), target       :: cache  ! Tree to walk trhough
      type(grCacheWalker), intent(inout) :: walker ! Walker structure
      !
      integer(ik) :: max_stack ! Maximum possible stack depth in walk-through
      integer(ik) :: alloc
      !
      max_stack = sum(cache%n_max+1)
      allocate (walker%n_stack(max_stack),walker%q_stack(max_stack),walker%quanta(size(cache%n_max)),stat=alloc)
      if (alloc/=0) stop 'cache_grx%walker_gr_initialize - allocation failed'
      walker%cache    => cache
      walker%node     => NULL() 
      walker%done     = .false.
      !$ call omp_init_lock(walker%traverse_lock)
    end subroutine walker_gr_initialize

    subroutine walker_gr_destroy(walker)
      type(grCacheWalker), intent(inout) :: walker ! Walker structure
      !
      if (allocated(walker%n_stack)) deallocate (walker%n_stack)
      if (allocated(walker%q_stack)) deallocate (walker%q_stack)
      if (allocated(walker%quanta))  deallocate (walker%quanta)
      nullify (walker%cache)
      !$ call omp_destroy_lock(walker%traverse_lock)
    end subroutine walker_gr_destroy

    function walker_gr_next(walker,n) result(success)
      type(grCacheWalker), intent(inout) :: walker  ! Walker structure
      integer(sik), intent(out)          :: n(:)    ! State vector, suitable for passing to cache_gr_locate()
      logical                            :: success ! .true. if n(:) contains the next node within the tree;
                                                    ! .false. if the walk is over
      !
      integer(ik) :: sp ! Stack pointer; syntax sugar for walker%p_stack
      integer(ik) :: pq ! Pointer to quanta counter in current model syntax sugar for walker%p_quanta
      integer(ik) :: iq, jq
      type(grCache), pointer :: dbg_cgr
      !
      if (size(n)/=size(walker%cache%n_max)) stop 'cache_grx%walker_gr_next - bad n(:) size'
      if (walker%done) then
        success = .false.
        return
      end if
      !$ call omp_set_lock(walker%traverse_lock)
      if (.not.associated(walker%node)) then
        !
        !  This is the first call to walker_gr_next(); start at the root of the tree
        !
        walker%node => walker%cache%root
        walker%p_quanta  = 1 ! We are working on the first mode
        walker%quanta(1) = 0 ! ,,, and it is not excited
        walker%p_stack   = 0 ! Nothing in the stack
      else
        sp = walker%p_stack
        pq = walker%p_quanta
        !
        !  We may have to back-track multiple times before we have an acceptable step
        !
        try_step: do
          if (associated(walker%node%left)) then
            !
            !  Try going to the left if we can; this increases the number of quanta in the leftmost mode
            !
            sp = sp + 1
            if (sp>size(walker%n_stack)) stop 'cache_grx%walker_gr_next - blown the stack (L)'
            walker%n_stack(sp)%p => walker%node
            walker%q_stack(sp)   = pq
            walker%node          => walker%node%left
            walker%quanta(pq)    = walker%quanta(pq) + 1
            exit try_step
          else
            !
            !  Ascend until we can step right
            !
            ascend: do while (.not.associated(walker%node%right)) 
              if (sp<=0) then
                !
                !  We've run our of nodes; return failure
                !
                walker%done = .true.
                success     = .false.
                !$ call omp_unset_lock(walker%traverse_lock)
                return
              end if
              walker%node       => walker%n_stack(sp)%p
              pq                = walker%q_stack(sp)
              walker%quanta(pq) = walker%quanta(pq) - 1
              sp                = sp - 1
            end do ascend
            !
            !  Step right now
            !
            walker%node       => walker%node%right
            pq                = pq + 1
            if (pq>size(walker%quanta)) stop 'cache_grx%walker_gr_next - blown modes array'
            walker%quanta(pq) = 0 
            !
            !  After a right step, we are not done yet: a pure right step gets us to 
            !  a clone of the the mode we stepped from (it has zero quanta in the
            !  mode we've just added). Therefore, we have to step left before we
            !  can accept the move.
            !
          end if
        end do try_step
        walker%p_stack  = sp
        walker%p_quanta = pq
      end if
      success = .true.
      !
      !  We now have a packed node index in walker%quanta(:pq). we have to expand it
      !  Since we delete zeros at the end of the compressed state vector, we can't
      !  just use unpack().
      !
      pq = walker%p_quanta
      jq = 0
      uncompress_state: do iq=1,size(n)
        if (jq>=pq .or. walker%cache%n_max(iq)==0) then
          n(iq) = 0
        else
          jq = jq + 1
          n(iq) = walker%quanta(jq)
        end if
      end do uncompress_state
      !$ call omp_unset_lock(walker%traverse_lock)
      !
      !  DEBUG: try to locate the current node by the state vector; we must get the 
      !  DEBUG: same node we are currently looking at!
      !
        dbg_cgr => cache_gr_locate(walker%cache,n)
        if (.not.associated(dbg_cgr,target=walker%node)) then
          write (out,"(('State vector: ',40i4))") n
          write (out,"('Expected node ',i16)") loc(walker%node)
          write (out,"(' Located node ',i16)") loc(dbg_cgr)
          stop 'cache_grx%walker_gr_next - bad node'
        end if
    end function walker_gr_next
    !
    !  Debug - dump tree pointers
    !
    recursive subroutine dump_tree_node(node,n_left,n_right)
      type(grCache), intent(in)   :: node
      integer(hik), intent(inout) :: n_left, n_right
      !
      ! write (out,"('Node ',i16,' left ',i16,' right ',i16)") loc(node), loc(node%left), loc(node%right)
      if (associated(node%left) ) then
        n_left = n_left + 1
        call dump_tree_node(node%left,n_left,n_right)
      end if
      if (associated(node%right)) then
        n_right = n_right + 1
        call dump_tree_node(node%right,n_left,n_right)
      end if
    end subroutine dump_tree_node
    !
    subroutine dump_tree(tree)
      type(grCacheWrapper), intent(in) :: tree ! Tree to dump
      integer(hik) :: n_left, n_right
      !
      ! write (out,"('=== begin tree connectivity dump ===')")
      n_left = 0 ; n_right = 0
      call dump_tree_node(tree%root,n_left,n_right)
      ! write (out,"('=== end tree connectivity dump ===')")
      write (out,"(' Number of left steps = ',i0)") n_left
      write (out,"('Number of right steps = ',i0)") n_right
    end subroutine dump_tree
    !
    !  Debug - complete tree walker
    !
    subroutine walker_complete(tree)
      type(grCacheWrapper), target :: tree ! Tree to walk
      !
      integer(hik)        :: node_count
      integer(sik)        :: n(size(tree%n_max))
      type(grCacheWalker) :: walker
      !
      call dump_tree(tree)
      call walker_gr_initialize(tree,walker)
      node_count = 0
      write (out,"()")
      scan_tree: do while(walker_gr_next(walker,n))
        node_count = node_count + 1
        write (out,"('Node ',i10,' is ',40i4)") node_count, n
      end do scan_tree
      call walker_gr_destroy(walker)
      write (out,"('Total number of nodes in the tree was ',i0/)") node_count
    end subroutine walker_complete
    !
  end module cache_grx
