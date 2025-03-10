!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_era5
! \label{get_era5}
!
!
! !REVISION HISTORY:
! 23 dec 2019: Sujay Kumar, initial code
!
! !INTERFACE:
subroutine get_era5(n, findex)
! !USES:
  use LDT_coreMod
  use LDT_timeMgrMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_metforcingMod
  use era5_forcingMod

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Opens, reads, and interpolates 1-hourly ERA5 forcing.
!
!  The ERA5 forcing data are organized into monthly files, where each
!  file contains 24 one-hourly records of forcing fields.  
!
!  In general, metforcing readers read the forcing data before the current
!  time, referred to as bookend1, and after the current time, referred to as
!  bookend2.  Then the readers temporally interpolate between bookend1 and
!  bookend2.  
!
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    forcing dataset index
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[LDT\_tick](\ref{LDT_tick}) \newline
!    call to advance or retract time
!  \item[era5files](\ref{era5files}) \newline
!    Puts together appropriate file name for 1 hour intervals
!  \item[read\_era5](\ref{read_era5}) \newline
!    call to read the ERA5 data and perform spatial interpolation
!  \end{description}
!EOP
  integer           :: order
  integer           :: ferror
  character(len=LDT_CONST_PATH_LEN) :: fname
  integer           :: c, r,kk,f,try
  integer           :: yr1, mo1, da1, hr1, mn1, ss1, doy1
  integer           :: yr2, mo2, da2, hr2, mn2, ss2, doy2
  real*8            :: time1, time2, timenow
  real*8            :: dtime1, dtime2
  real              :: gmt1, gmt2
  real              :: ts1, ts2

  integer           :: hr_int1, hr_int2
  integer           :: movetime  ! Flag to move bookend2 files to bookend1

! _________________________________________________________

  if( LDT_rc%nts(n).gt.3600 ) then   ! > 1-hr timestep
     write(LDT_logunit,*) '[ERR] When running LDT with ERA5, the clock '
     write(LDT_logunit,*) '[ERR] should run with a timestep less than or '
     write(LDT_logunit,*) '[ERR] equal to one hour.'
     call LDT_endrun()
  endif

  era5_struc(n)%findtime1 = 0
  era5_struc(n)%findtime2 = 0
  movetime = 0

!=== Determine Required ERA5 Data Times (The previous hour and the future hour)
  yr1=LDT_rc%yr
  mo1=LDT_rc%mo
  da1=LDT_rc%da
  hr1=LDT_rc%hr
  mn1=LDT_rc%mn
  ss1=0
  ts1=0
  call LDT_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
 
  if(LDT_rc%ts.gt.3600) then 
     write(LDT_logunit,*) '[ERR] The model timestep is > forcing data timestep'
     write(LDT_logunit,*) '[ERR] LDT does not support this mode currently'
     write(LDT_logunit,*) '[ERR] Program stopping ...'
     call LDT_endrun()
  endif

  if(mod(nint(LDT_rc%ts),3600).eq.0) then 
     if(timenow.ge.era5_struc(n)%era5time2) then 
        yr1 = LDT_rc%yr
        mo1=LDT_rc%mo
        da1=LDT_rc%da
        hr1=LDT_rc%hr
        mn1=0
        ss1=0
        ts1=-60*60
        call LDT_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        
        yr2=LDT_rc%yr    !next hour
        mo2=LDT_rc%mo
        da2=LDT_rc%da
        hr2=LDT_rc%hr
        mn2=0
        ss2=0
        ts2=0
        call LDT_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        movetime = 1
        era5_struc(n)%findtime2 = 1
     endif
  else
     if(timenow.ge.era5_struc(n)%era5time2) then 
        yr1 = LDT_rc%yr
        mo1=LDT_rc%mo
        da1=LDT_rc%da
        hr1=LDT_rc%hr
        mn1=0
        ss1=0
        ts1=0
        call LDT_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

        yr2=LDT_rc%yr    !next hour
        mo2=LDT_rc%mo
        da2=LDT_rc%da
        hr2=LDT_rc%hr
        mn2=0
        ss2=0
        ts2=60*60
        call LDT_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

        movetime = 1
        era5_struc(n)%findtime2 = 1
     endif
  endif

  if(LDT_rc%tscount(n).eq.1 .or.LDT_rc%rstflag(n).eq.1  ) then  
     era5_struc(n)%findtime1=1
     era5_struc(n)%findtime2=1
     movetime=0
     LDT_rc%rstflag(n) = 0
  endif
  
  if(movetime.eq.1) then
     era5_struc(n)%era5time1=era5_struc(n)%era5time2
     do f=1,LDT_rc%met_nf(findex)
        do c=1,LDT_rc%ngrid(n)
           era5_struc(n)%metdata1(f,c)=era5_struc(n)%metdata2(f,c)
        enddo
     enddo
  endif    !end of movetime=1
  
  if(era5_struc(n)%findtime1.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferror=0
     try=0  
     ts1=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1
        kk = 1
        order = 1
        call era5files(n,kk,findex,era5_struc(n)%era5dir, yr1, mo1, da1, &
             fname)
        write(unit=LDT_logunit,fmt=*)'[INFO] getting file.. ',trim(fname)
        call read_era5(n, kk,order, yr1,mo1, da1,hr1,&
             findex, fname, ferror)

        if(ferror.ge.1) era5_struc(n)%era5time1=time1
        call LDT_tick(dtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        if(try.gt.11)then
           write(LDT_logunit,*)'[ERR] ERA5 data gap exceeds 10 days on file 1'
           call LDT_endrun()
        endif
     enddo
!=== end of data search
  endif   

  if(era5_struc(n)%findtime2.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferror=0
     try=0  
     ts2=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1

     !- Obtaining ERA5 File:
        kk=1
        order = 2
        call era5files(n,kk,findex,era5_struc(n)%era5dir, yr2, mo2, da2, &
             fname)
        write(unit=LDT_logunit,fmt=*)'[INFO] getting file.. ',trim(fname)
        call read_era5(n, kk,order, yr2,mo2,da2,hr2, &
             findex, fname, ferror)

        if(ferror.ge.1) then
           era5_struc(n)%era5time2=time2
        endif
        call LDT_tick(dtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        if(try.gt.11)then
           write(LDT_logunit,*)'[ERR] ERA5 data gap exceeds 10 days on file 2'
           call LDT_endrun()
        endif
     enddo
  endif 
end subroutine get_era5


!BOP
! !ROUTINE: era5files
! \label{era5files}
!
! !INTERFACE:
subroutine era5files(n, kk, findex, era5dir, yr, mo, da, fname)

! !USES:
  use LDT_coreMod
  use LDT_logMod
  use LDT_timeMgrMod

  implicit none
! !ARGUMENTS:
  integer                       :: n 
  integer                       :: kk
  integer                       :: findex
  character(len=*), intent(in)  :: era5dir
  integer, intent(in)           :: yr,mo,da
  character(len=*), intent(out) :: fname

! !DESCRIPTION:
!   This subroutine puts together ERA5 file names for
!   daily netcdf files
!
!  The arguments are:
!  \begin{description}
!  \item[era5dir]
!    Name of the ERA5 directory
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[fname]
!   name of the timestamped ERA5 file
!  \end{description}
!
!EOP

  character*4  :: cyear
  character*2  :: cmonth
  character*8  :: cdate
  integer      :: hr, mn, ss
  real*8       :: time
  integer      :: doy
  real         :: gmt

  hr = 0 
  mn = 0 
  ss = 0 

  write(unit=cyear, fmt='(i4.4)') yr
  write(unit=cmonth,fmt='(i2.2)') mo
  
  fname = trim(era5dir)//'/FORCING_'//trim(cyear)//trim(cmonth)//'.nc'
end subroutine era5files

