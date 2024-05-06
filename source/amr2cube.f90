program amr2cube
  !--------------------------------------------------------------------------
  ! Ce programme calcule le cube cartesien pour les
  ! variables hydro d'une simulation RAMSES. 
  ! Version F90 par R. Teyssier le 01/04/01.
  !--------------------------------------------------------------------------
  use ramses_info
  use hydro_reader
  implicit none
  character(128)::type
  integer::n,i,j,k,ncoarse
  integer::nvar,ncpuh,ngrid_current
  integer::nx,ny,nz,ilevel,idim,jdim,kdim,icell
  integer::ilevel1,ngrid1
  integer::nlevelmaxs,nlevel,iout
  integer::ind,ipos,ngridh,ilevela,ilevelh
  integer::nhx,nhy,ihx,ihy,ivar1,ivar2
  real::smallr,smallc,gammah
  real::boxlen,boxlen2
  real::hexp,t2,aexp2,hexp2
  real::scale_l,scale_d,scale_t
  real::omega_m2,omega_l2,omega_k2,omega_b2

  integer::nx_sample=0,ny_sample=0,nz_sample=0
  integer::ngrid,imin,imax,jmin,jmax,kmin,kmax
  integer::ncpu2,npart2,ndim2,nlevelmax2,nstep_coarse2
  integer::nx2,ny2,nz2,ngridmax2,ndimh,nlevelmaxh
  integer::nx_full,ny_full,nz_full,lmin
  integer::ix,iy,iz,ixp1,iyp1,izp1,ndom,impi,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8)::order_min
  !real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax
  real(KIND=8)::ddx,ddy,ddz,dex,dey,dez,xx,yy,zz
!  real(KIND=8),dimension(:,:),allocatable::x
  !==========================================
  real(KIND=8),dimension(:,:),allocatable::pos
  !==========================================
  real(KIND=4),dimension(:,:,:),allocatable::toto
  real(KIND=8),dimension(:)  ,allocatable::rho
  real(KIND=4)::x_off=0,y_off=0,z_off=0,dx_header=0
  integer,dimension(:)  ,allocatable::isp
  character(LEN=5)::ncharcpu  
  character(LEN=1000)::nomfich,outfich,filetype='bin'
  character(LEN=1000)::repository
  logical::ok,ok_part,ok_cell
  character(LEN=1)::proj='z'
  logical::do_readgrav=.false.
  logical::do_readkrome=.false.
  logical::readRT=.true.
  type level
     integer::ilevel
     integer::ngrid
     real(KIND=4),dimension(:,:,:),pointer::cube
     integer::imin
     integer::imax
     integer::jmin
     integer::jmax
     integer::kmin
     integer::kmax
  end type level

  type(level),dimension(1:100)::grid


  call read_params
  call read_info(repository)
  call prepareEmissivity(type)
  call init_reader(repository, readRT, do_readgrav, do_readkrome)


  if(ndim==2)then
     write(*,*)'Output file contains 2D data'
     write(*,*)'Aborting'
     stop
  endif

  !-----------------------
  ! Map parameters
  !-----------------------
  xxmin=xmin ; xxmax=xmax
  yymin=ymin ; yymax=ymax
  zzmin=zmin ; zzmax=zmax

  !-----------------------------
  ! Compute hierarchy
  !-----------------------------
  do ilevel=1,lmax
     nx_full=2**ilevel
     ny_full=2**ilevel
     nz_full=2**ilevel
     imin=int(xxmin*dble(nx_full))+1
     imax=int(xxmax*dble(nx_full))+1
     jmin=int(yymin*dble(ny_full))+1
     jmax=int(yymax*dble(ny_full))+1
     kmin=int(zzmin*dble(nz_full))+1
     kmax=int(zzmax*dble(nz_full))+1
     allocate(grid(ilevel)%cube(imin:imax,jmin:jmax,kmin:kmax))
     grid(ilevel)%cube=0.0
     grid(ilevel)%imin=imin
     grid(ilevel)%imax=imax
     grid(ilevel)%jmin=jmin
     grid(ilevel)%jmax=jmax
     grid(ilevel)%kmin=kmin
     grid(ilevel)%kmax=kmax
  end do

  !-----------------------------------------------
  ! Compute projected variables
  !----------------------------------------------

  ! Loop over processor files
  print*,'Number of cpus to read = ',ncpu_read
  do k=1,ncpu_read
     call read_cpu(k,repository,ncpu_read)

     ! Loop over levels
     do ilevel=1,lmax
        write (*, "(A, f5.2, A, A)", advance='no') &           ! Progress bar
             ' Progress ',dble(k) / ncpu_read * 100,' % ',char(13)
        call read_level(ilevel)
        ! Geometry
        nx_full=2**ilevel
        ny_full=2**ilevel
        nz_full=2**ilevel

        ! Allocate work arrays
        grid(ilevel)%ngrid=ngrida
        if(ngrida>0)then
           allocate(rho(1:ngrida))
        endif

        ! Compute map
        if(ngrida>0)then

           ! Loop over cells
           do ind=1,twotondim

              ! Extract variable
              if(index(type, 'cpu') .eq. 1) then !------------------------
                 rho = icpu
              else if(index(type, 'ref') .eq. 1) then !-------------------
                 rho = ilevel
              else
                 call emissivity(var, rho, ind, dx)
              endif

              ! Store data cube
              do i=1,ngrida
                 ok_cell= .not.ref(i,ind)
                 !ok_cell=.true.
                 if(ok_cell)then
                    ix=int(xp(i,ind,1)*dble(nx_full))+1
                    iy=int(xp(i,ind,2)*dble(ny_full))+1
                    iz=int(xp(i,ind,3)*dble(nz_full))+1
                    if(    ix>=grid(ilevel)%imin.and.&
                         & iy>=grid(ilevel)%jmin.and.&
                         & iz>=grid(ilevel)%kmin.and.&
                         & ix<=grid(ilevel)%imax.and.&
                         & iy<=grid(ilevel)%jmax.and.&
                         & iz<=grid(ilevel)%kmax)then
                       grid(ilevel)%cube(ix,iy,iz)=rho(i)
                    endif
                 end if
              end do

           end do
           ! End loop over cell

           deallocate(rho)
        endif

     end do
     ! End loop over levels

     call close_cpu()

  end do
  ! End loop over cpus

  nx_full=2**lmax
  ny_full=2**lmax
  nz_full=2**lmax
  imin=int(xxmin*dble(nx_full))+1
  imax=int(xxmax*dble(nx_full))
  jmin=int(yymin*dble(ny_full))+1
  jmax=int(yymax*dble(ny_full))
  kmin=int(zzmin*dble(nz_full))+1
  kmax=int(zzmax*dble(nz_full))

  write(*,'(" Zoom in",3(1x,I5,"-",I5))')imin,imax,jmin,jmax,kmin,kmax
  write(*,'(" res=",3(1x,I5))')imax-imin+1,jmax-jmin+1,kmax-kmin+1
  write(*,*)'Please wait...'
  
  do ix=imin,imax
     xmin=((dble(ix)-0.5d0)/2**lmax)
     do iy=jmin,jmax
        ymin=((dble(iy)-0.5d0)/2**lmax)
        do iz=kmin,kmax
           zmin=((dble(iz)-0.5d0)/2**lmax)
           do ilevel=max(levelmin,lmin),lmax-1
              ndom=2**ilevel
              i=int(xmin*ndom)+1
              j=int(ymin*ndom)+1
              k=int(zmin*ndom)+1
              grid(lmax)%cube(ix,iy,iz)=grid(lmax)%cube(ix,iy,iz)+ &
                   & grid(ilevel)%cube(i,j,k)
           end do
        end do
     end do
  end do
  write(*,*)'Min value:',minval(grid(lmax)%cube(imin:imax,jmin:jmax,kmin:kmax))
  write(*,*)'Max value:',maxval(grid(lmax)%cube(imin:imax,jmin:jmax,kmin:kmax))
  write(*,*)'Norm:     ',sum   (dble(grid(lmax)%cube(imin:imax,jmin:jmax,kmin:kmax))) &
       & /(imax-imin+1)/(jmax-jmin+1)/(kmax-kmin+1)

  ! Output file
  if(TRIM(filetype).eq.'bin')then
     nomfich=TRIM(outfich)
     write(*,*)'Writing file '//TRIM(nomfich)
     open(unit=20,file=nomfich,form='unformatted')
     if(nx_sample==0)then
        write(20)imax-imin+1,jmax-jmin+1,kmax-kmin+1
        allocate(toto(imax-imin+1,jmax-jmin+1,kmax-kmin+1))
        toto=grid(lmax)%cube(imin:imax,jmin:jmax,kmin:kmax)
        print*,'min and max values = ',minval(toto),maxval(toto)
        !toto=log10(toto)
        write(20)toto
     else
        if(ny_sample==0)ny_sample=nx_sample
        if(nz_sample==0)nz_sample=ny_sample
        write(20)nx_sample,ny_sample,nz_sample
        allocate(toto(1:nx_sample,1:ny_sample,1:nz_sample))
        do i=1,nx_sample
           xx=(dble(i)-0.5)/dble(nx_sample)*dble(imax-imin+1)
           ix=int(xx)
           ddx=xx-ix
           dex=1d0-ddx
           ix=ix+imin
           ix=min(ix,imax)
           ixp1=min(ix+1,imax)
           do j=1,ny_sample
              yy=(dble(j)-0.5)/dble(ny_sample)*dble(jmax-jmin+1)
              iy=int(yy)
              ddy=yy-iy
              dey=1d0-ddy
              iy=iy+jmin
              iy=min(iy,jmax)
              iyp1=min(iy+1,jmax)
              do k=1,nz_sample
                 zz=(dble(k)-0.5)/dble(nz_sample)*dble(kmax-kmin+1)
                 iz=int(zz)
                 ddz=zz-iz
                 dez=1d0-ddz
                 iz=iz+kmin
                 iz=min(iz,kmax)
                 izp1=min(iz+1,kmax)
                 toto(i,j,k)=dex*dey*dez*log10(grid(lmax)%cube(ix  ,iy  ,iz  ))+ &
                      &      ddx*dey*dez*log10(grid(lmax)%cube(ixp1,iy  ,iz  ))+ &
                      &      dex*ddy*dez*log10(grid(lmax)%cube(ix  ,iyp1,iz  ))+ &
                      &      dex*dey*ddz*log10(grid(lmax)%cube(ix  ,iy  ,izp1))+ &
                      &      ddx*ddy*dez*log10(grid(lmax)%cube(ixp1,iyp1,iz  ))+ &
                      &      dex*ddy*ddz*log10(grid(lmax)%cube(ix  ,iyp1,izp1))+ &
                      &      ddx*dey*ddz*log10(grid(lmax)%cube(ixp1,iy  ,izp1))+ &
                      &      ddx*ddy*ddz*log10(grid(lmax)%cube(ixp1,iyp1,izp1))
              end do
           end do
        end do
        toto=10.**(toto)
        write(20)toto
     endif
     write(20)xxmin,xxmax
     write(20)yymin,yymax
     write(20)zzmin,zzmax
     close(20)
  endif
  
  if(TRIM(filetype).eq.'grafic')then
     nomfich=TRIM(outfich)
     write(*,*)'Writing file '//TRIM(nomfich)
     open(unit=20,file=nomfich,form='unformatted')
     dummy=0.0
     if (dx_header==0) then   
         write(20)int(imax-imin+1),int(jmax-jmin+1),int(kmax-kmin+1),1./(imax-imin),x_off,y_off,z_off,dummy,dummy,dummy,dummy
      else
         write(20)int(imax-imin+1),int(jmax-jmin+1),int(kmax-kmin+1),dx_header,x_off,y_off,z_off,dummy,dummy,dummy,dummy
         print*, 'header :', dx_header,x_off,y_off,z_off
      endif
     do iz=kmin,kmax
        write(20)((real(grid(lmax)%cube(ix,iy,iz)),ix=imin,imax),iy=jmin,jmax)
     end do
     close(20)
  endif
  
  if(TRIM(filetype).eq.'vtk')then
     nomfich=TRIM(outfich)
     write(*,*)'Writing file '//TRIM(nomfich)
     open(unit=20,file=nomfich,form='formatted')
     write(20,'("# vtk DataFile Version 2.0")')
     write(20,'("RAMSES data using vtk file format")')
     write(20,'("ASCII")')
     write(20,'("DATASET STRUCTURED_POINTS")')
     write(20,'("DIMENSIONS ",3(I3,1x))')imax-imin+1,jmax-jmin+1,kmax-kmin+1
     write(20,'("ORIGIN 0.0 0.0 0.0")')
     write(20,'("SPACINGS 1.0 1.0 1.0")')
     write(20,'(" ")')
     write(20,'("POINT_DATA ",I8)')(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)
     write(20,'("SCALARS values float")')
     write(20,'("LOOKUP_TABLE default")')
     do k=kmin,kmax
        do j=jmin,jmax
           do i=imin,imax
              write(20,'(F9.3)')real((grid(lmax)%cube(i,j,k)))
           end do
        end do
     end do
     close(20)
  endif
  
contains
  
  subroutine read_params
    
    implicit none
    
    integer       :: i,n
    integer       :: iargc
    character(len=4)   :: opt
    character(len=1000) :: arg
    !LOGICAL       :: bad, ok
    
    n = iargc()
    if (n < 4) then
       print *, 'usage: amr2cube  -inp  input_dir'
       print *, '                 -out  output_file'
       print *, '                 [-xmi xmin] '
       print *, '                 [-xma xmax] '
       print *, '                 [-ymi ymin] '
       print *, '                 [-yma ymax] '
       print *, '                 [-zmi zmin] '
       print *, '                 [-zma zmax] '
       print *, '                 [-lma lmax] '
       print *, '                 [-typ type] '
       print *, '                 [-fil filetype] '
       print *, 'ex: amr2cube -inp output_00001 -out cube.dat'// &
            &   ' -typ 1 -xmi 0.1 -xma 0.7 -lma 12'
       !print *, ' '
       !print *, ' type :-1 = cpu number'
       !print *, ' type : 0 = ref. level (default)'
       !print *, ' type : 1-9 = variable number'
       stop
    end if
    
    do i = 1,n,2
       call getarg(i,opt)
       if (i == n) then
          print '("option ",a2," has no argument")', opt
          stop 2
       end if
       call getarg(i+1,arg)
       select case (opt)
       case ('-inp')
          repository = trim(arg)
       case ('-out')
          outfich = trim(arg)
       case ('-fil')
          filetype = trim(arg)
       case ('-xmi')
          read (arg,*) xmin
       case ('-xma')
          read (arg,*) xmax
       case ('-ymi')
          read (arg,*) ymin
       case ('-yma')
          read (arg,*) ymax
       case ('-zmi')
          read (arg,*) zmin
       case ('-zma')
          read (arg,*) zmax
       case ('-lma')
          read (arg,*) lmax
       case ('-nx')
          read (arg,*) nx_sample
       case ('-ny')
          read (arg,*) ny_sample
       case ('-nz')
          read (arg,*) nz_sample
       case ('-typ')
          type = trim(arg)
       case ('-xof')
          read (arg,*) x_off
       case ('-yof')
          read (arg,*) y_off
       case ('-zof')
          read (arg,*) z_off
       case ('-dx')
          read (arg,*) dx_header
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do
    
    return
    
  end subroutine read_params
  
end program amr2cube

!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character*5::nchar

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif


end subroutine title
!================================================================
!================================================================
!================================================================
!================================================================
