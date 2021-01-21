!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_gdasnc
! \label{read_gdasnc}
!
! !REVISION HISTORY:
!  14 Dec 2000: Urszula Jambor; Rewrote geteta.f to use GDAS forcing in GLDAS
!  15 Mar 2001: Jon Gottschalck; Added additional parameters and octets in 
!               which to search in GRIB files
!  01 Jun 2001: Urszula Jambor; Added option to get forcing from different 
!               files (F00 instantaneous and F06 six hour means)
!  29 Jan 2003: Urszula Jambor; Rewrote code, uses GETGB call to replace
!               ungribgdas.  Interpolation now occurs in interp_gdasnc.  
!               Using GETGB avoids problems with the Oct2002 GDAS 
!               grid update.
!  12 Nov 2003: Matt Rodell; Check to make sure input file exists before
!		opening and thereby creating a new, empty file.
!  14 Nov 2003: Matt Rodell; Ensure lugb varies in call to baopen
!  05 Feb 2004: James Geiger; Added GrADS-DODS Server functionality
!  29 Apr 2010: Sujay Kumar: Fixed the mixing of instantaneous and time 
!               averaged fields.
!  25 Jan 2012: Sujay Kumar; Switched to the use of grib-api library
!  11 Jan 2021: KR Arsenault; added new GDAS NetCDF-4 formatted file support
!
! !INTERFACE:
subroutine read_gdasnc( order, n, findex, &
     name00, name03, name06, F06flag, ferror,try )
! !USES:  
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_timeMgrMod,     only : LIS_date2time
  use LIS_metforcingMod,  only : LIS_forc
  use gdasnc_forcingMod,  only : gdas_struc
  use LIS_logMod,         only : LIS_logunit, LIS_endrun
  use LIS_surfaceModelDataMod

  implicit none
! !ARGUMENTS:
  integer, intent(in)          :: order    
  integer, intent(in)          :: n
  integer, intent(in)          :: findex
  character(len=*), intent(in) :: name00
  character(len=*), intent(in) :: name03
  character(len=*), intent(in) :: name06
  logical, intent(in)          :: F06flag
  integer, intent(out)         :: ferror 
  integer, intent(inout)       :: try
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  GDAS forecast datasets, transforms into 9 LIS forcing 
!  parameters and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read for the previous 
!    3hr bookend, order=2, read for the next 3 hr bookend)
!  \item[n]
!    index of the nest
!  \item[name00]
!    name of the instantaneous analysis file
!  \item[name03]
!    name of the 3 hour GDAS forecast file
!  \item[name06]
!    name of the 6 hour GDAS forecast file
!  \item[F06flag]
!    flag to indicate if 6hr forecast data is required for this interval
!  \item[ferror]
!    flag to indicate success of the call (=0 indicates success)
!  \item[try]
!    index of the tries (in case of missing data)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_gdas](\ref{interp_gdasnc}) \newline
!    spatially interpolates a GDAS variable
!  \end{description}
!EOP
!==== Local Variables=======================
  
  character(120) :: filename
  integer :: iv, c,r,t
  integer :: ferror1, ferror2, ferror3
  integer :: ngdas
  real    :: glbdata_i(10,LIS_rc%ngrid(n))
  real    :: glbdata_a(10,LIS_rc%ngrid(n))
  real    :: glbdata_a_f06(10,LIS_rc%ngrid(n))
  logical :: dataStrucflag


!=== End Variable Definition =======================

  glbdata_i = LIS_rc%udef
  glbdata_a = LIS_rc%udef
  glbdata_a_f06 = LIS_rc%udef
  ngdas = (gdas_struc(n)%ncold*gdas_struc(n)%nrold)
  dataStrucflag = .false.

!--------------------------------------------------------------------------
! Check model timestep and output interval and stop if the timestep
! and/or output interval will cause incorrect output to be generated.
!--------------------------------------------------------------------------
  if(LIS_rc%ts .ge. 10800 .AND. LIS_rc%ts .lt. 21600) then
     write(LIS_logunit,*) "[ERR] Model timestep is greater than or equal"
     write(LIS_logunit,*) "[ERR]   to 3 hours, which will cause errors in the "
     write(LIS_logunit,*) "[ERR]   GDAS reader. Change the model timestep to a"
     write(LIS_logunit,*) "[ERR]   value less than 3 hours ('1hr' is suggested.)"
     call LIS_endrun()
 endif
  if(LIS_sfmodel_struc(n)%outInterval .ge. 21600 .AND. LIS_rc%ts .eq. 21600) then
     write(LIS_logunit,*) '[ERR] Model timestep is 6hr and output interval is greater than'
     write(LIS_logunit,*) '[ERR]   or equal to 6hr. This setup can cause issues in the reader'
     write(LIS_logunit,*) '[ERR]   where the output is not the true 6hr average and/or the '
     write(LIS_logunit,*) '[ERR]   data being written is shifted by one timestep.'
     write(LIS_logunit,*) '[ERR] It is suggested that the model timestep be changed to 1hr or less.'
     call LIS_endrun()
  endif

!--------------------------------------------------------------------------
! If there's a problem then ferror is set to zero
! read instantaneous fields
!--------------------------------------------------------------------------
  iv = 0

!--------------------------------------------------------------------------
! Set up to open file and retrieve specified field 
!--------------------------------------------------------------------------
  filename = name00
  if( gdas_struc(n)%dstrucchange1 .and. &
      gdas_struc(n)%gdastime1 .ge. gdas_struc(n)%datastructime1 ) then
    dataStrucflag = .true.  ! HKB Use special routine for f00 files following 2019 Jun 12 12Z GDAS upgrades
  endif
  call retrieve_gdasnc_variables(n, findex, filename, dataStrucflag, glbdata_i, ferror1)
  dataStrucflag = .false.   ! Reset flag since f03 and f06 files are not affected by 2019 Jun 12 upgrade

!--------------------------------------------------------------------------
! Read 3hr forecast for time averaged fields
!--------------------------------------------------------------------------
  filename = name03
  call retrieve_gdasnc_variables(n, findex, filename, dataStrucflag, glbdata_a, ferror2)

!--------------------------------------------------------------------------
! Read 6hr forecast for time averaged fields, if required. 
!--------------------------------------------------------------------------
  if(F06flag) then 
     filename = name06
     call retrieve_gdasnc_variables(n, findex, filename, dataStrucflag, &
                                    glbdata_a_f06, ferror3)
  end if
  
  ferror = 1
  if(ferror1.eq.0.or.ferror2.eq.0) then 
     ferror = 0
  endif
  if(F06flag) then
     if(ferror.eq.0.or.ferror3.eq.0) then 
        ferror = 0
     endif
  endif
!--------------------------------------------------------------------------
! Place the interpolated data into the LIS arrays
!--------------------------------------------------------------------------
  do iv=1,9
     do t=1,LIS_rc%ngrid(n)
        if(F06flag) then 
           if(iv.eq.3.or.iv.eq.4.or.iv.eq.8.or.iv.eq.9) then ! these are time avgd fields
              if(order.eq.1) then 
                 gdas_struc(n)%metdata1(iv,t) = 2*glbdata_a_F06(iv,t)-glbdata_a(iv,t)
              else
                 gdas_struc(n)%metdata2(iv,t) = 2*glbdata_a_F06(iv,t)-glbdata_a(iv,t)
              endif
           else ! these are instantaneous
              if(order.eq.1) then 
                 gdas_struc(n)%metdata1(iv,t) = glbdata_i(iv,t)
              else
                 gdas_struc(n)%metdata2(iv,t) = glbdata_i(iv,t)
              endif
           endif
        else
           if(iv.eq.3.or.iv.eq.4.or.iv.eq.8.or.iv.eq.9) then ! these are time avgd fields
              if(order.eq.1) then 
                 gdas_struc(n)%metdata1(iv,t) = glbdata_a(iv,t)                         
              else
                 gdas_struc(n)%metdata2(iv,t) = glbdata_a(iv,t)                         
              endif
           else
              if(order.eq.1) then 
                 gdas_struc(n)%metdata1(iv,t) = glbdata_i(iv,t)                         
              else
                 gdas_struc(n)%metdata2(iv,t) = glbdata_i(iv,t)                         
              endif
           endif
        endif
     enddo
  enddo

end subroutine read_gdasnc

!BOP
! 
! !ROUTINE: retrieve_gdasnc_variables
! \label{retrieve_gdasnc_variables}
! 
! !INTERFACE: 
subroutine retrieve_gdasnc_variables(n, findex, filename, dataStrucflag,&
                                     glbdata, errorcode)
! !USES: 
  use LIS_coreMod,       only : LIS_rc, LIS_domain
  use LIS_logMod,        only : LIS_logunit, LIS_endrun,& 
                           LIS_verify, LIS_warning
  use gdasnc_forcingMod, only : gdas_struc

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer              :: n      ! Nest
  integer              :: findex
  character(len=*)     :: filename
  logical              :: dataStrucflag
  real                 :: glbdata(10,LIS_rc%ngrid(n))
  integer              :: errorcode
! 
! !DESCRIPTION: 
!  This subroutine retrieves GDAS-Netcdf4 forcing variables, and 
!   interpolates them to the LIS grid. 
! 
!EOP
  integer           :: ngdas
  integer           :: iv,ivmax
  integer           :: c,r,t
  real              :: missingValue 
  real, allocatable :: f(:)       ! Input 1-D forcing field
  real, allocatable :: f2d(:,:)   ! Input 2-D forcing field
  real, dimension(LIS_rc%lnc(n), LIS_rc%lnr(n)) :: varfield

! Netcdf-4 parameters and fields:
  integer   :: ftn, varId, ncId, nrId
  integer   :: numcols, numrows
  character(14), dimension(9), parameter :: gdasnc_fv = (/  &
       'TMP_2m_inst   ',  &
       'SPFH_2m_inst  ',  &
       'DSWRF_sfc_tavg',  &
       'DLWRF_sfc_tavg',  &
       'UGRD_10m_inst ',  &
       'VGRD_10m_inst ',  &
       'PRES_sfc_inst ',  &
       'PRATE_sfc_tavg',  &
       'CPRAT_sfc_tavg'   /)

! - End of Netcdf parameter / var declarations:

  logical :: file_exists
  integer :: var_index
  logical :: pcp_flag
  logical*1, allocatable :: lb(:)

! _________________________________________

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  if(dataStrucflag) then
    ! HKB...All instantaneous fields in f00 files
!    pds16 = (/010,010,010,010,010,010,010,010,010 /)
    ivmax = 7
  else
    ! index 10 indicates instantaneous, 003 indicates time average
!    pds16 = (/010,010,003,003,010,010,010,003,003 /) 
    ivmax = 9
  endif

  ngdas = (gdas_struc(n)%ncold*gdas_struc(n)%nrold)

  varfield = 0 
  errorcode = 1

  ! See if file exists:
  inquire (file=filename, exist=file_exists)
  if (file_exists) then      

     ! Open file, if exists:
     call LIS_verify(nf90_open(path=trim(filename), mode=NF90_NOWRITE, &
                     ncId=ftn), 'nf90_open failed in read_gdasnc')

     call LIS_verify( nf90_get_att(ftn, NF90_GLOBAL, 'missing_value', missingValue), &
                    'Error in nf90_get_att for missing_value in read_gdasnc')

     ! Read in lat / lon dimensions to check against grids within current date/time:
     call LIS_verify( nf90_inq_dimid(ftn,"lon",ncId), &
                    'Error in nf90_inq_dimid for lon in read_gdasnc')
     call LIS_verify( nf90_inq_dimid(ftn,"lat",nrId), &
                    'Error in nf90_inq_dimid for lat in read_gdasnc')
 
     call LIS_verify( nf90_inquire_dimension(ftn, ncId, len=numcols), &
                    'Error in nf90_inquire_dimension for numcols in read_gdasnc')
     call LIS_verify( nf90_inquire_dimension(ftn, nrId, len=numrows), &
                    'Error in nf90_inquire_dimension for numrows in read_gdasnc')

     ! Double check that input file rows and cols match the reader values:
     if( ngdas /= (numcols*numrows) ) then
       write(LIS_logunit,*) &
           '[ERR] Number of values does not match those in: ',trim(filename)
       write(LIS_logunit,*) &
           '[ERR] Number of points in GDAS-NC reader: ',ngdas
       write(LIS_logunit,*) &
           '[ERR] Number of points in input GDAS-NC file: ',(numcols*numrows)
       call LIS_endrun 
     endif

     allocate(lb(gdas_struc(n)%ncold*gdas_struc(n)%nrold))
     allocate(f(gdas_struc(n)%ncold*gdas_struc(n)%nrold))
     allocate(f2d(gdas_struc(n)%ncold,gdas_struc(n)%nrold)) 

     ! Loop over GDAS-netcdf variables:
     var_index = 0
     do iv = 1, gdas_struc(n)%nmif

        call LIS_verify( nf90_inq_varid(ftn, trim(gdasnc_fv(iv)), varid), &
                'nf90_inq_varid failed for var read in read_gdasnc')

        ! Reading in 2D GDAS forcing field:
        f2d = LIS_rc%udef
        call LIS_verify( nf90_get_var(ftn, varid, f2d, &
                         start=(/1,1/), &
                         count=(/gdas_struc(n)%ncold,gdas_struc(n)%nrold/)),&
                'nf90_get_var failed for var read in read_gdasnc')
! KRA TO BE UPDATED TO SUPPORT PARALLEL READING IN OF FILES ...

        var_index = iv 
        lb = .false.
        do t=1,ngdas
           if( f(t).ne.missingValue ) then
             lb(t) = .true. 
           endif
        enddo
           
        pcp_flag = .false. 
        if( var_index.eq.8.or.var_index.eq.9 ) then 
           pcp_flag = .true. 
        endif

        ! Convert 2d var array to 1d var array for interpolation step:
        t = 0 
        f = LIS_rc%udef
        do r = 1, gdas_struc(n)%nrold
           do c = 1, gdas_struc(n)%ncold
              t = t + 1
              f(t) = f2d(c,r)
           enddo 
        enddo
           
        ! Pass in 1D array of each met variable to interpolate
        !  to LIS local parallel subdomain grid:
        call interp_gdasnc(n, findex, pcp_flag, ngdas, f, &
                    lb, LIS_rc%gridDesc(n,:), &
                    LIS_rc%lnc(n), LIS_rc%lnr(n), varfield )
           
        ! Convert interpolated 2D grid back to 1D:
        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              if( LIS_domain(n)%gindex(c,r).ne.-1 ) then 
                  glbdata(var_index,LIS_domain(n)%gindex(c,r)) = &
                         varfield(c,r)
              endif
           enddo
        enddo

     enddo  ! END VARIABLE LOOP 

     ! Close netcdf file:
     call LIS_verify(nf90_close(ftn))

     deallocate(lb)
     deallocate(f)     
         
  else
     write(LIS_logunit,*) &
          '[ERR] Could not find file: ',trim(filename)
     errorcode = 0
  endif
! End of retrieving GDAS NetCDF-4 fields
#endif     

end subroutine retrieve_gdasnc_variables


!BOP
! !ROUTINE: interp_gdasnc
! \label{interp_gdasnc}
!
! !INTERFACE:
subroutine interp_gdasnc(n, findex, pcp_flag, ngdas, &
                         f, lb, lis_gds, nc, nr, varfield)
! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_domain
  use LIS_logMod,        only : LIS_logunit, LIS_endrun
  use gdasnc_forcingMod, only : gdas_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
  logical, intent(in) :: pcp_flag
  integer, intent(in) :: ngdas
  real, intent(out)   :: f(ngdas)
  logical*1           :: lb(ngdas)
  real                :: lis_gds(50)
  integer, intent(in) :: nc
  integer, intent(in) :: nr
  real, intent(out)   :: varfield(nc,nr)
!
! !DESCRIPTION:
!   This subroutine interpolates a given GDAS field 
!   to the LIS grid. 
!
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[ngdas]
!  number of elements in the input grid
! \item[f]
!  input data array to be interpolated
! \item[lb]
!  input bitmap
! \item[lis\_gds]
!  array description of the LIS grid
! \item[nc]
!  number of columns (in the east-west dimension) in the LIS grid
! \item[nr]
!  number of rows (in the north-south dimension) in the LIS grid
! \item[varfield]
!  output interpolated field
!  \end{description} 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!   \item[upscaleByAveraging](\ref{upscaleByAveraging}) \newline
!     upscales finer resolution forcing data to coarser resolution running
!     domain by averaging
! \end{description}
!EOP
  integer :: iret
  integer :: count1,i,j,mo
  real, dimension(nc*nr) :: lis1d
  logical*1 :: lo(nc*nr)

!=== End variable declarations
!-----------------------------------------------------------------------
! Setting interpolation options (ip=0,bilinear)
! (km=1, one parameter, ibi=1,use undefined bitmap
! (needed for soil moisture and temperature only)
! Use budget bilinear (ip=3) for precip forcing fields
!-----------------------------------------------------------------------
  mo = nc*nr
!-----------------------------------------------------------------------
! Initialize output bitmap. Important for soil moisture and temp.
!-----------------------------------------------------------------------
  lo = .true.
!-----------------------------------------------------------------------  
! Interpolate to LIS grid
!-----------------------------------------------------------------------  
  if ( gdas_struc(n)%met_interp == "bilinear" ) then
     call bilinear_interp(lis_gds,lb,f,lo,lis1d,gdas_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          gdas_struc(n)%w111,gdas_struc(n)%w121,&
          gdas_struc(n)%w211,gdas_struc(n)%w221,&
          gdas_struc(n)%n111,gdas_struc(n)%n121,&
          gdas_struc(n)%n211,gdas_struc(n)%n221,LIS_rc%udef, iret)
  elseif ( gdas_struc(n)%met_interp == "budget-bilinear" ) then
     if (pcp_flag) then 
        call conserv_interp(lis_gds,lb,f,lo,lis1d,gdas_struc(n)%mi,mo, & 
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             gdas_struc(n)%w112,gdas_struc(n)%w122,&
             gdas_struc(n)%w212,gdas_struc(n)%w222,&
             gdas_struc(n)%n112,gdas_struc(n)%n122,&
             gdas_struc(n)%n212,gdas_struc(n)%n222,LIS_rc%udef,iret)
     else 
        call bilinear_interp(lis_gds,lb,f,lo,lis1d,gdas_struc(n)%mi,mo,&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             gdas_struc(n)%w111,gdas_struc(n)%w121,&
             gdas_struc(n)%w211,gdas_struc(n)%w221,&
             gdas_struc(n)%n111,gdas_struc(n)%n121,&
             gdas_struc(n)%n211,gdas_struc(n)%n221,LIS_rc%udef,iret)
     endif

  elseif ( gdas_struc(n)%met_interp == "average" ) then
    call upscaleByAveraging(gdas_struc(n)%mi, mo, LIS_rc%udef, &
         gdas_struc(n)%n111, lb, f, lo, lis1d)
  else
     write(LIS_logunit,*) '[ERR] The specified spatial interpolation option '
     write(LIS_logunit,*) '[ERR] is not supported for GDAS.'
     write(LIS_logunit,*) '[ERR] LIS is stopping.'
     call LIS_endrun()
  endif

!-----------------------------------------------------------------------    
! Create 2D array for main program. Also define a "soil" mask
! due to different geography between GDAS & LDAS. For LDAS land 
! points not included in GDAS geography dataset only.
!-----------------------------------------------------------------------    
  count1 = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = lis1d(i+count1)
     enddo
     count1 = count1 + nc
  enddo

end subroutine interp_gdasnc

