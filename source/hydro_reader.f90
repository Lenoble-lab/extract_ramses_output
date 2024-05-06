MODULE hydro_reader
  !--------------------------------------------------------------------------
  ! Read amr and hydro structure into variables
  ! By Joki Rosdahl, 21/10/2013
  !--------------------------------------------------------------------------
  use ramses_info
  implicit none
  integer::nboundary
  integer::ngrida,ncpu_read,lmax,twotondim,icpu,nvarh,nvarRT

  real(KIND=8),dimension(:,:),allocatable::xg
  real(KIND=8),dimension(:,:,:),allocatable::xp ! Cell center positions
  real(KIND=8),dimension(:,:,:),allocatable::var
  real(KIND=4),dimension(:),allocatable::var_sp

  logical,dimension(:,:),allocatable::ref
  integer,dimension(:,:),allocatable::son,ngridfile,ngridlevel,ngridbound

  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,dmax,dummy
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
  real(KIND=8)::dx
  real(KIND=8),dimension(1:8,1:3)::xc
  real(KIND=8),dimension(1:3)::xbound=(/0d0,0d0,0d0/)
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list

  integer,parameter::qdp=kind(1.0_16) ! real*16
  real(kind=qdp),dimension(:),allocatable::bound_key
  
  logical::rtOK=.false.,gravOK=.false.,krome=.false.

CONTAINS

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE init_reader(snsnap,do_readRT,do_readgrav,do_readkrome)

! Prepare reading of amr structure
!-------------------------------------------------------------------------
  character(LEN=128)::nomfich,snsnap
  logical,optional::do_readRT, do_readgrav, do_readkrome
  logical::ok,readRT=.false.,readgrav=.false.
  real::boxlen
  integer::nx,ny,nz,ilevel,impi,i,j
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,order_min
  integer::bit_length,lmin,ndom,maxdom,ngrid_current
  integer::imin,imax,jmin,jmax,kmin,kmax
!-------------------------------------------------------------------------
  nomfich=TRIM(snsnap)//'/hydro_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     !stop
  endif
  ! Try RT variables:
  if(present(do_readRT)) readRT=do_readRT
  if(readRT) then
     nomfich=TRIM(snsnap)//'/rt_'//TRIM(nchar)//'.out00001'
     inquire(file=nomfich, exist=ok) ! verify input file
     if ( .not. ok ) then
        print *,TRIM(nomfich)//' not found.'
        rtOk=.false.
     else 
        rtOK=.true.
     endif
  endif
  ! Read gravity variables:
  if(present(do_readgrav)) readgrav=do_readgrav
  if(readgrav) then
     nomfich=TRIM(snsnap)//'/grav_'//TRIM(nchar)//'.out00001'
     inquire(file=nomfich, exist=ok) ! verify input file 
     if ( .not. ok ) then
        print *,TRIM(nomfich)//' not found.'
        !stop
     endif
     gravOK=.true.
  endif
  
  !read krome variable if krome
  if(present(do_readkrome)) krome=do_readkrome


  nomfich=TRIM(snsnap)//'/amr_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif

  nomfich=TRIM(snsnap)//'/amr_'//TRIM(nchar)//'.out00001'
  open(unit=10,file=nomfich,status='old',form='unformatted')
  read(10)ncpu
  read(10)ndim
  read(10)nx,ny,nz
  read(10)nlevelmax
  read(10)ngridmax
  read(10)nboundary
  read(10)ngrid_current
  read(10)boxlen
  close(10)
  twotondim=2**ndim
  xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)
  if(allocated(ngridfile)) deallocate(ngridfile, ngridlevel, cpu_list)
  allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
  allocate(ngridlevel(1:ncpu,1:nlevelmax))
  if(allocated(ngridbound)) deallocate(ngridbound)
  if(nboundary>0)allocate(ngridbound(1:nboundary,1:nlevelmax))

  !if(ndim==2)then
  !   write(*,*)'Output file contains 2D data'
  !   write(*,*)'Aborting'
  !   stop
  !endif

  allocate(cpu_list(1:ncpu))
  if(TRIM(ordering).eq.'hilbert')then
     if(allocated(cpu_read)) deallocate(cpu_read)
     allocate(cpu_read(1:ncpu))
     cpu_read=.false.
  endif

  !:::Leo: from RASCAS (devs by Maxime Trebitsch and LMD)
  if(allocated(bound_key)) deallocate(bound_key)
  allocate(bound_key(0:ncpu))
  call read_hilbert_keys_raw(snsnap,ncpu,bound_key)
  !:::oeL
  
  !-----------------------
  ! Map parameters
  !-----------------------
  if(lmax==0)then
     lmax=nlevelmax
  endif
  !write(*,*)'time, level and resolution=',t,lmax,2**lmax
  !write(*,*)'Working resolution =',2**lmax
  xxmin=xmin ; xxmax=xmax
  yymin=ymin ; yymax=ymax
  zzmin=zmin ; zzmax=zmax

  if(TRIM(ordering).eq.'hilbert')then

     dmax=max(xxmax-xxmin,yymax-yymin,zzmax-zzmin)
     do ilevel=1,lmax
        dx=0.5d0**ilevel
        if(dx.lt.dmax*4)exit
     end do
     lmin=ilevel
     bit_length=lmin-1
     maxdom=2**bit_length
     imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
     if(bit_length>0)then
        imin=int(xxmin*dble(maxdom))
        imax=imin+1
        jmin=int(yymin*dble(maxdom))
        jmax=jmin+1
        kmin=int(zzmin*dble(maxdom))
        kmax=kmin+1
     endif

     dkey=(dble(2**(nlevelmax+1)/dble(maxdom)))**ndim
     ndom=1
     if(bit_length>0)ndom=8
     idom(1)=imin; idom(2)=imax
     idom(3)=imin; idom(4)=imax
     idom(5)=imin; idom(6)=imax
     idom(7)=imin; idom(8)=imax
     jdom(1)=jmin; jdom(2)=jmin
     jdom(3)=jmax; jdom(4)=jmax
     jdom(5)=jmin; jdom(6)=jmin
     jdom(7)=jmax; jdom(8)=jmax
     kdom(1)=kmin; kdom(2)=kmin
     kdom(3)=kmin; kdom(4)=kmin
     kdom(5)=kmax; kdom(6)=kmax
     kdom(7)=kmax; kdom(8)=kmax
     
     do i=1,ndom
        if(bit_length>0)then
!!$           call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
           call hilbert3d_1pt(idom(i),jdom(i),kdom(i),order_min,bit_length)
        else
           order_min=0.0d0
        endif
        bounding_min(i)=(dble(order_min))*dkey
        bounding_max(i)=(dble(order_min)+1.0D0)*dkey
     end do
     
     cpu_min=0; cpu_max=0
     do impi=1,ncpu
        do i=1,ndom
           if (   bound_key(impi-1).le.bounding_min(i).and.&
                & bound_key(impi  ).gt.bounding_min(i))then
              cpu_min(i)=impi
           endif
           if (   bound_key(impi-1).lt.bounding_max(i).and.&
                & bound_key(impi  ).ge.bounding_max(i))then
              cpu_max(i)=impi
           endif
        end do
     end do
     
     ncpu_read=0
     do i=1,ndom
        do j=cpu_min(i),cpu_max(i)
           if(.not. cpu_read(j))then
              ncpu_read=ncpu_read+1
              cpu_list(ncpu_read)=j
              cpu_read(j)=.true.
           endif
        enddo
     enddo
  else
     ncpu_read=ncpu
     do j=1,ncpu
        cpu_list(j)=j
     end do
  end  if

END SUBROUTINE init_reader

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE read_cpu(kcpu,repo,ncpu_read)

! Prepare reading of amr structure
!-------------------------------------------------------------------------
  integer::kcpu,i,ncpu_read
  character(LEN=5)::ncharcpu
  !character(LEN=3)::percDone
  character(LEN=128)::nomfich,repo
!-------------------------------------------------------------------------
  icpu=cpu_list(kcpu)
  call title(icpu,ncharcpu)

  ! Open AMR file and skip header
  nomfich=TRIM(repo)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
  open(unit=10,file=nomfich,status='old',form='unformatted')
  !write(*,*)'Processing file '//TRIM(nomfich)
  !if(mod(kcpu*10,ncpu_read) .lt. 10)  then
  !   write(percDone,'(I2)') (kcpu*10)/ncpu_read
  !   write(*,"(A)",advance="no") percDone
  !endif
  do i=1,21
     read(10)
  end do
  ! Read grid numbers
  read(10)ngridlevel
  ngridfile(1:ncpu,1:nlevelmax)=ngridlevel
  read(10)
  if(nboundary>0)then
     do i=1,2
        read(10)
     end do
     read(10)ngridbound
     ngridfile(ncpu+1:ncpu+nboundary,1:nlevelmax)=ngridbound
  endif
  read(10)
  ! ROM: comment the single follwing line for old stuff
  read(10)
  if(TRIM(ordering).eq.'bisection')then
     do i=1,5
        read(10)
     end do
  else
     read(10)
  endif
  read(10)
  read(10)
  read(10)
  
  ! Open HYDRO file and skip header
  nomfich=TRIM(repo)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
  open(unit=11,file=nomfich,status='old',form='unformatted')
  read(11)
  read(11)nvarh
  read(11)
  read(11)
  read(11)
  read(11)

  !add 2 variables if krome activate (for mu and Tgas)
  if (krome) then
      nvarh = nvarh + 2
  endif

  nvarRT=0
  if (rtOk) then
     ! Open RT file and skip header
     nomfich=TRIM(repo)//'/rt_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=12,file=nomfich,status='old',form='unformatted')
     read(12)
     read(12)nvarRT
     read(12)
     read(12)
     read(12)
     read(12)
  endif

  if (gravOK) then
     ! Open GRAV file and skip header
     nomfich=TRIM(repo)//'/grav_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=13,file=nomfich,status='old',form='unformatted')
     read(13)
     read(13)
     read(13)
     read(13)
  endif

END SUBROUTINE read_cpu
  
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE read_level(ilevel)

! Read level from already opened repository
!-------------------------------------------------------------------------
  integer::ilevel,ind,idim,ivar
  integer::ix,iy,iz,i,j
!-------------------------------------------------------------------------
  ! Geometry
  dx=0.5**ilevel
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
     
  ! Allocate work arrays
  if(allocated(xg)) then 
     deallocate(xg,son,var,xp,ref)
  endif
  if(allocated(var_sp)) deallocate(var_sp)
  ngrida=ngridfile(icpu,ilevel)
!!$  grid(ilevel)%ngrid=ngrida
  if(ngrida>0)then
     allocate(xg(1:ngrida,1:ndim))
     allocate(son(1:ngrida,1:twotondim))
     if(gravOK)then
        allocate(var(1:ngrida,1:twotondim,1:nvarh+nvarRT+1+ndim))
     else
        allocate(var(1:ngrida,1:twotondim,1:nvarh+nvarRT))
     endif
     if((rtOK .and. rtPrecision.eq.4) .or. Uprecision.eq.4) allocate(var_sp(1:ngrida))
     allocate(xp(1:ngrida,1:twotondim,1:ndim))
     allocate(ref(1:ngrida,1:twotondim))
     ref=.false.
  endif
  
  ! Loop over domains
  do j=1,nboundary+ncpu
     
     ! Read AMR data
     if(ngridfile(j,ilevel)>0)then
        read(10) ! Skip grid index
        read(10) ! Skip next index
        read(10) ! Skip prev index
        ! Read grid center
        do idim=1,ndim
           if(j.eq.icpu)then
              read(10)xg(:,idim)
           else
              read(10)
           endif
        end do
        read(10) ! Skip father index
        do ind=1,2*ndim
           read(10) ! Skip nbor index
        end do
        ! Read son index
        do ind=1,twotondim
           if(j.eq.icpu)then
              read(10)son(:,ind)
           else
              read(10)
           end if
        end do
        ! Skip cpu map
        do ind=1,twotondim
           read(10)
        end do
        ! Skip refinement map
        do ind=1,twotondim
           read(10)
        end do
     endif
     
     ! Read HYDRO data
     read(11)
     read(11)
     if(rtOk)read(12)
     if(rtOk)read(12)
     if(gravOK)read(13)
     if(gravOK)read(13)
     if(ngridfile(j,ilevel)>0)then
        ! Read hydro variables
        do ind=1,twotondim
           do ivar=1,nvarh
              if(j.eq.icpu)then
                 if(Uprecision.eq.4) then
                    read(11) var_sp(:)
                    var(:,ind,ivar) = var_sp(:)
                 else
                    read(11)var(:,ind,ivar)
                 endif
              else
                 read(11)
              end if
           end do
           do ivar=1,nvarRT
              if(j.eq.icpu)then
                 if(rtPrecision.eq.4) then
                    read(12) var_sp(:)
                    var(:,ind,nvarh+ivar) = var_sp(:)
                 else
                    read(12)var(:,ind,nvarh+ivar)
                 endif
              else
                 read(12)
              end if
           end do
           ! Read grav potential and force
           if(gravOK)then
              ! Read potential
              if(j.eq.icpu)then
                 if(Uprecision.eq.4) then
                    read(13) var_sp(:)
                    var(:,ind,nvarh+nvarRT+1) = var_sp(:)
                 else
                    read(13)var(:,ind,nvarh+nvarRT+1)
                 endif
              else
                 read(13)
              end if
              ! Read force
              do ivar=1,ndim
                 if(j.eq.icpu)then
                    if(Uprecision.eq.4) then
                       read(13) var_sp(:)
                       var(:,ind,nvarh+nvarRT+1+ivar) = var_sp(:)
                    else
                       read(13)var(:,ind,nvarh+nvarRT+1+ivar)
                    endif
                 else
                    read(13)
                 endif
              enddo
           endif
        end do
     end if
  end do

  ! Compute cell center positions
  if(ngrida>0)then
     ! Loop over cells
     do ind=1,twotondim
        ! Compute cell center
        do i=1,ngrida
           xp(i,ind,1)=(xg(i,1)+xc(ind,1)-xbound(1))
           xp(i,ind,2)=(xg(i,2)+xc(ind,2)-xbound(2))
           if(ndim.gt.2) xp(i,ind,3)=(xg(i,3)+xc(ind,3)-xbound(3))
        end do      
        ! Check if cell is refined
        do i=1,ngrida
           ref(i,ind)=son(i,ind)>0.and.ilevel<lmax
        end do
     end do     
  endif


END SUBROUTINE read_level
  

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE close_cpu()

! 
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  close(10)
  close(11)
  if(rtOk) close(12)
  if(gravOK) close(13)
END SUBROUTINE close_cpu  
!================================================================
!================================================================
!================================================================
!================================================================
subroutine read_hilbert_keys_raw(snsnap,ncpu,bound_key)
  ! read hilbert keys in one amr file
  ! from RASCAS
  
  implicit none
  character(128),intent(in)                 :: snsnap
  integer(kind=4),intent(in)                :: ncpu
  real(qdp),dimension(0:ncpu),intent(inout) :: bound_key
  real(kind=8),dimension(0:ncpu)            :: bound_key_dp    
  logical(kind=4)                           :: ok
  character(128)                            :: nomfich
  character(128)                            :: orderingtype
  integer(kind=4)                           :: i, ios
  integer(kind=4),parameter                 :: param_unit = 13

  write(nomfich,'(a,a,a,a)') trim(snsnap), '/amr_', trim(nchar), '.out00001'
  inquire(file=nomfich, exist=ok)
  if(.not. ok)then
     write(*,*)'File '//TRIM(nomfich)//' not found'    
     stop
  end if
  open(unit=param_unit,file=nomfich,form='unformatted',status='old',action='read',iostat=ios)
  do i=1,24 ! Assume that there is no "simple boundary"
     read(param_unit,iostat=ios)
  end do
  read(param_unit,iostat=ios) orderingtype

  if (trim(orderingtype) .ne. 'bisection') then
     read(param_unit,iostat=ios) bound_key
     if(ios/=0) then ! read in dp
        close(param_unit)
        open(unit=param_unit,file=nomfich,form='unformatted',status='old',action='read',iostat=ios)
        do i=1,24 ! Assume that there is no "simple boundary"
           read(param_unit,iostat=ios)
        end do
        read(param_unit,iostat=ios) orderingtype
        read(param_unit,iostat=ios) bound_key_dp
        if(ios/=0) then
           call read_hilbert_keys(snsnap,ncpu,bound_key_dp)
           print*,'Hilbert keys read in the info file'
        else
           print*,'Hilbert keys read in double precision'
        end if
        bound_key = real(bound_key_dp, kind=qdp)
     else
        print*,'Hilbert keys read in quad precision'
     end if
  end if

  close (param_unit)

  return

end subroutine read_hilbert_keys_raw
!================================================================
!================================================================
!================================================================
!================================================================
subroutine read_hilbert_keys(snsnap,ncpu,bound_key)
  ! read the hilbert keys in the info file
  ! from RASCAS
  
  implicit none
  character(128),intent(in)                    :: snsnap
  integer(kind=4),intent(in)                   :: ncpu
  real(kind=8),dimension(0:ncpu),intent(inout) :: bound_key
  logical(kind=4)            :: not_ok
  character(512)             :: filename
  character(512)             :: line,name,value,orderingtype
  integer(kind=4)            :: i, impi
  integer(kind=4),parameter  :: param_unit = 13

  not_ok = .true.
  write(filename,'(a,a,a,a)') trim(snsnap),'/info_',trim(nchar),'.txt'
  open(unit=param_unit,file=filename,status='old',form='formatted')
  do
     read(param_unit,'(a)',end=2) line
     i = scan(line,'=')
     if (i==0 .or. line(1:1)=='#') cycle
     name=trim(adjustl(line(:i-1)))
     value=trim(adjustl(line(i+1:)))
     ! check for a comment at end of line !
     i = scan(value,'!')
     if (i /= 0) value = trim(adjustl(value(:i-1)))
     if (trim(name) .eq. 'ordering type') then
        read(value,*) orderingtype
        if (orderingtype=='hilbert') then
           ! start reading keys
           ! skip one line
           read(param_unit,*)
           do impi=1,ncpu
              read(param_unit,'(I8,1X,E23.15,1X,E23.15)')i,bound_key(impi-1),bound_key(impi)
           end do
        end if
     end if
  end do
2 close (param_unit)
  
  return
end subroutine read_hilbert_keys


END MODULE hydro_reader


!================================================================
!================================================================
!================================================================
!================================================================
subroutine hilbert3d(x,y,z,order,bit_length,npoint)
  implicit none

  integer     ,INTENT(IN)                     ::bit_length,npoint
  integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
  real(kind=8),INTENT(OUT),dimension(1:npoint)::order

  logical,dimension(0:3*bit_length-1)::i_bit_mask
  logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
  integer,dimension(0:7,0:1,0:11)::state_diagram
  integer::i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit

  if(bit_length>bit_size(bit_length))then
     write(*,*)'Maximum bit length=',bit_size(bit_length)
     write(*,*)'stop in hilbert3d'
     stop
  endif

  state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
                            &   0, 1, 3, 2, 7, 6, 4, 5,&
                            &   2, 6, 0, 7, 8, 8, 0, 7,&
                            &   0, 7, 1, 6, 3, 4, 2, 5,&
                            &   0, 9,10, 9, 1, 1,11,11,&
                            &   0, 3, 7, 4, 1, 2, 6, 5,&
                            &   6, 0, 6,11, 9, 0, 9, 8,&
                            &   2, 3, 1, 0, 5, 4, 6, 7,&
                            &  11,11, 0, 7, 5, 9, 0, 7,&
                            &   4, 3, 5, 2, 7, 0, 6, 1,&
                            &   4, 4, 8, 8, 0, 6,10, 6,&
                            &   6, 5, 1, 2, 7, 4, 0, 3,&
                            &   5, 7, 5, 3, 1, 1,11,11,&
                            &   4, 7, 3, 0, 5, 6, 2, 1,&
                            &   6, 1, 6,10, 9, 4, 9,10,&
                            &   6, 7, 5, 4, 1, 0, 2, 3,&
                            &  10, 3, 1, 1,10, 3, 5, 9,&
                            &   2, 5, 3, 4, 1, 6, 0, 7,&
                            &   4, 4, 8, 8, 2, 7, 2, 3,&
                            &   2, 1, 5, 6, 3, 0, 4, 7,&
                            &   7, 2,11, 2, 7, 5, 8, 5,&
                            &   4, 5, 7, 6, 3, 2, 0, 1,&
                            &  10, 3, 2, 6,10, 3, 4, 4,&
                            &   6, 1, 7, 0, 5, 2, 4, 3 /), &
                            & (/8 ,2, 12 /) )

  do ip=1,npoint

     ! convert to binary
     do i=0,bit_length-1
        x_bit_mask(i)=btest(x(ip),i)
        y_bit_mask(i)=btest(y(ip),i)
        z_bit_mask(i)=btest(z(ip),i)
     enddo

     ! interleave bits
     do i=0,bit_length-1
        i_bit_mask(3*i+2)=x_bit_mask(i)
        i_bit_mask(3*i+1)=y_bit_mask(i)
        i_bit_mask(3*i  )=z_bit_mask(i)
     end do

     ! build Hilbert ordering using state diagram
     cstate=0
     do i=bit_length-1,0,-1
        b2=0 ; if(i_bit_mask(3*i+2))b2=1
        b1=0 ; if(i_bit_mask(3*i+1))b1=1
        b0=0 ; if(i_bit_mask(3*i  ))b0=1
        sdigit=b2*4+b1*2+b0
        nstate=state_diagram(sdigit,0,cstate)
        hdigit=state_diagram(sdigit,1,cstate)
        i_bit_mask(3*i+2)=btest(hdigit,2)
        i_bit_mask(3*i+1)=btest(hdigit,1)
        i_bit_mask(3*i  )=btest(hdigit,0)
        cstate=nstate
     enddo

     ! save Hilbert key as double precision real
     order(ip)=0.
     do i=0,3*bit_length-1
        b0=0 ; if(i_bit_mask(i))b0=1
        order(ip)=order(ip)+dble(b0)*dble(2)**i
     end do

  end do

end subroutine hilbert3d

subroutine hilbert3d_1pt(x,y,z,order,bit_length)
  implicit none

  integer     ,INTENT(IN)  ::bit_length
  integer     ,INTENT(IN)  ::x,y,z
  real(kind=8),INTENT(OUT) ::order

  logical,dimension(0:3*bit_length-1)::i_bit_mask
  logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
  integer,dimension(0:7,0:1,0:11)::state_diagram
  integer::i,cstate,nstate,b0,b1,b2,sdigit,hdigit

  if(bit_length>bit_size(bit_length))then
     write(*,*)'Maximum bit length=',bit_size(bit_length)
     write(*,*)'stop in hilbert3d'
     stop
  endif

  state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
                            &   0, 1, 3, 2, 7, 6, 4, 5,&
                            &   2, 6, 0, 7, 8, 8, 0, 7,&
                            &   0, 7, 1, 6, 3, 4, 2, 5,&
                            &   0, 9,10, 9, 1, 1,11,11,&
                            &   0, 3, 7, 4, 1, 2, 6, 5,&
                            &   6, 0, 6,11, 9, 0, 9, 8,&
                            &   2, 3, 1, 0, 5, 4, 6, 7,&
                            &  11,11, 0, 7, 5, 9, 0, 7,&
                            &   4, 3, 5, 2, 7, 0, 6, 1,&
                            &   4, 4, 8, 8, 0, 6,10, 6,&
                            &   6, 5, 1, 2, 7, 4, 0, 3,&
                            &   5, 7, 5, 3, 1, 1,11,11,&
                            &   4, 7, 3, 0, 5, 6, 2, 1,&
                            &   6, 1, 6,10, 9, 4, 9,10,&
                            &   6, 7, 5, 4, 1, 0, 2, 3,&
                            &  10, 3, 1, 1,10, 3, 5, 9,&
                            &   2, 5, 3, 4, 1, 6, 0, 7,&
                            &   4, 4, 8, 8, 2, 7, 2, 3,&
                            &   2, 1, 5, 6, 3, 0, 4, 7,&
                            &   7, 2,11, 2, 7, 5, 8, 5,&
                            &   4, 5, 7, 6, 3, 2, 0, 1,&
                            &  10, 3, 2, 6,10, 3, 4, 4,&
                            &   6, 1, 7, 0, 5, 2, 4, 3 /), &
                            & (/8 ,2, 12 /) )

  ! convert to binary
  do i=0,bit_length-1
     x_bit_mask(i)=btest(x,i)
     y_bit_mask(i)=btest(y,i)
     z_bit_mask(i)=btest(z,i)
  enddo

  ! interleave bits
  do i=0,bit_length-1
     i_bit_mask(3*i+2)=x_bit_mask(i)
     i_bit_mask(3*i+1)=y_bit_mask(i)
     i_bit_mask(3*i  )=z_bit_mask(i)
  end do
  
  ! build Hilbert ordering using state diagram
  cstate=0
  do i=bit_length-1,0,-1
     b2=0 ; if(i_bit_mask(3*i+2))b2=1
     b1=0 ; if(i_bit_mask(3*i+1))b1=1
     b0=0 ; if(i_bit_mask(3*i  ))b0=1
     sdigit=b2*4+b1*2+b0
     nstate=state_diagram(sdigit,0,cstate)
     hdigit=state_diagram(sdigit,1,cstate)
     i_bit_mask(3*i+2)=btest(hdigit,2)
     i_bit_mask(3*i+1)=btest(hdigit,1)
     i_bit_mask(3*i  )=btest(hdigit,0)
     cstate=nstate
  enddo

  ! save Hilbert key as double precision real
  order=0.
  do i=0,3*bit_length-1
     b0=0 ; if(i_bit_mask(i))b0=1
     order=order+dble(b0)*dble(2)**i
  end do

end subroutine hilbert3d_1pt

!================================================================
!================================================================
!================================================================
!================================================================
