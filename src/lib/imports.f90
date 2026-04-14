module imports
!*Brief Description:* This module contains all the subroutines required to
!import fields, previous model results, etc.
!*LICENSE:*
!
!
!
!*Full Description:*
!
  !
  use arrays
  use diagnostics
  use geometry
  use indices
  use other_consts
  use ventilation
  
  implicit none

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private
  public import_capillary
  public import_ventilation
  public import_perfusion
  public import_terminal
  
contains
  !
  !###########################################################################################
  !
  !>*import_capillary:* This subroutine reads in the results for the micro-circulatory
  ! components (capillary bed within 'units') of a perfusion model .
  subroutine import_capillary(micro_unit_file)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_IMPORT_CAPILLARY" :: IMPORT_CAPILLARY
    
    character(len=MAX_FILENAME_LEN),intent(in) :: micro_unit_file
    !local variables
    integer :: count_units,ios,ne,nunit
    real(dp) :: Pin,Pout,Qtot,temp(7),TOTAL_CAP_VOL,TOTAL_SHEET_SA,TT_TOTAL
    
    character(len=60) :: sub_name
    
    sub_name = 'import_capillary'
    call enter_exit(sub_name,1)
    
    open(10, file = micro_unit_file, status='old')
    
    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.
    ios = 0
    count_units = 0
    do while (ios == 0)
       read(10, '(I6,X,5(F9.2,X),2(F8.5,X),F10.2,X,F8.4,X,F10.4,X,F10.3,X,F8.4,X,F9.4,X)', iostat=ios)  &
            ne,temp(1:3),Pin,Pout,temp(4),Qtot,temp(5),TOTAL_CAP_VOL,TOTAL_SHEET_SA,TT_TOTAL,temp(6:7)
       if(ios == 0)then
          ! record the unit values for mean pressure, transit time, surface area.
          ! ne is the 'linker' element, so nunit is for its parent element
          nunit = int(elem_field(ne_unit,elem_cnct(-1,1,ne))) 
          unit_field(nu_blood_press,nunit) = (Pin+Pout)/2.0_dp
          unit_field(nu_tt,nunit) = TT_TOTAL
          unit_field(nu_sa,nunit) = TOTAL_SHEET_SA
          count_units = count_units + 1
       endif
    enddo

    if(count_units.ne.num_units)then
       write(*,'(''WARNING: the number of capillary units ('',i6,'') does not match the geometric model ('',i6,'')'')') &
            count_units,num_units
    endif

    close(10)
    
   call enter_exit(sub_name,2)
   
 end subroutine import_capillary
!
!##############################################################################
!
!>*import_ventilation:* This subroutine reads in the results of a ventilation model that
! has been saved in an exelem format as a single flow field (elements listed with
! ventilation as field values).
 subroutine import_ventilation(FLOWFILE)
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_IMPORT_VENTILATION" :: IMPORT_VENTILATION

   character(len=MAX_FILENAME_LEN),intent(in) :: FLOWFILE
   !local variables
   integer :: ierror,ne,nunit
   character(LEN=132) :: ctemp1,exfile
   real(dp) :: flow,flow_unit,maxflow

   character(len=60) :: sub_name

   sub_name = 'import_ventilation'
   call enter_exit(sub_name,1)

   print *, 'Reading in ventilation results'
   call import_exelemfield(FLOWFILE,ne_Vdot)
   do nunit = 1,num_units
     ne = units(nunit)
     if(elem_field(ne_Vdot,ne).lt.0.0_dp) elem_field(ne_Vdot,ne) = zero_tol
     unit_field(nu_Vdot0,nunit) = elem_field(ne_Vdot,ne)
   enddo

!!! sum the fields up the tree
   call sum_elem_field_from_periphery(ne_Vdot) !sum the air flows recursively UP the tree
   maxflow = elem_field(ne_Vdot,1)

   call enter_exit(sub_name,2)
 end subroutine import_ventilation

!
!###########################################################################################
!
!>*import_perfusion:* This subroutine reads in the results of a ventilation model that
! has been saved in an exelem format as a single flow field (elements listed with
! ventilation as field values).
 subroutine import_perfusion(FLOWFILE)
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_IMPORT_PERFUSION" :: IMPORT_PERFUSION

   character(len=MAX_FILENAME_LEN),intent(in) :: FLOWFILE
   !local variables
   integer :: ierror,ne,nunit
   character(LEN=132) :: ctemp1,exfile
   real(dp) :: flow,flow_unit,maxflow

   character(len=60) :: sub_name

   sub_name = 'import_perfusion'
   call enter_exit(sub_name,1)

   print *, 'Reading in perfusion results'
   call import_exelemfield(FLOWFILE,ne_Qdot)
   do nunit = 1,num_units
     ne = units(nunit)
     if(elem_field(ne_Qdot,ne).lt.0.0_dp) elem_field(ne_Qdot,ne) = zero_tol
     unit_field(nu_perf,nunit) = elem_field(ne_Qdot,ne)
   enddo

!!! sum the fields up the tree
   call sum_elem_field_from_periphery(ne_Qdot) !sum the air flows recursively UP the tree
   maxflow = elem_field(ne_Qdot,1)

   call enter_exit(sub_name,2)
 end subroutine import_perfusion

!
!##############################################################################
!
!>*import_exelemfield:* This subroutine reads in the content of an exelem field file (1 field)
 subroutine import_exelemfield(FLOWFILE,field_no)
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_IMPORT_EXELEMFIELD" :: IMPORT_EXELEMFIELD

   character(len=MAX_FILENAME_LEN),intent(in) :: FLOWFILE
   integer, intent(in) :: field_no
   !local variables
   integer :: ierror,ne,nunit
   character(LEN=132) :: ctemp1,exfile
   real(dp) :: flow,flow_unit,maxflow

   character(len=60) :: sub_name

   sub_name = 'import_exelemfield'
   call enter_exit(sub_name,1)

   open(10, file=FLOWFILE, status='old')
   ne = 0
   read_elem_flow : do !define a do loop name
     !.......read element flow
     read(unit=10, fmt="(a)", iostat=ierror) ctemp1
     if(index(ctemp1, "Values:")> 0) then
       ne = ne+1
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       flow = get_final_real(ctemp1)
       if(flow.lt.0.0_dp) flow = zero_tol
         elem_field(field_no,ne) = flow! read it in
       end if

!     if (ne .eq. ne+1)
!       print *, "ERROR: Refined element number does not match with mesh."
!     call exit(0)
       if(ne.ge.num_elems) exit read_elem_flow
     end do read_elem_flow

   close(10)

    call enter_exit(sub_name,2)
 end subroutine import_exelemfield

!
!##############################################################################
!
!>*import_terminal:* This subroutine reads in the content of an exnode field file
 subroutine import_terminal(EXFILE)
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_IMPORT_TERMINAL" :: IMPORT_TERMINAL

   character(len=MAX_FILENAME_LEN),intent(in) :: EXFILE
   !local variables
   integer :: count_units,field_label(10),i,ibeg,iend,ierror,ne,nunit,n_fields
   real(dp) :: rtemp
   character(LEN=132) :: ctemp,label
   character(len=300) :: readfile

   character(len=60) :: sub_name

   sub_name = 'import_terminal'
   call enter_exit(sub_name,1)

   if(index(EXFILE, ".exnode")> 0) then !full filename is given
      readfile = EXFILE
   else ! need to append the correct filename extension
      readfile = trim(EXFILE)//'.exnode'
   endif
   
   open(10, file=readfile, status='old')

   n_fields = 0
   ierror = 0

   read_field_labels : do
      read(10, fmt="(a)", iostat=ierror) ctemp
      if(index(ctemp, "Node:")> 0) exit read_field_labels
      if(index(ctemp, ") ")> 0) then
         ibeg = index(ctemp, ") ")+1 ! beginning of label
         iend = index(ctemp, ",")-1  ! end of label
         label = adjustl(ctemp(ibeg:iend))
         n_fields = n_fields + 1
         field_label(n_fields) = 0
         
         select case(label)
         case('terminal_element')
         case('pleural_pressure')
            field_label(n_fields) = nu_pe
         case('tidal_volume')
            field_label(n_fields) = nu_vt
         case('max_Pe')
            field_label(n_fields) = nu_Pe_max
         case('min_Pe')
            field_label(n_fields) = nu_Pe_min
         end select
      endif
         
   end do read_field_labels

   nunit = 0
   ierror = 0
   
   do while (ierror == 0)
      if(index(ctemp, "Node:")> 0) then
         
         do i = 1,3  ! read coordinates; not used 
            read(10, *, iostat=ierror) ctemp
         enddo
         read(10, *, iostat=ierror) ctemp ! read element: not used but could be

         nunit = nunit + 1
         
         do i = 3,n_fields
            read(10, *, iostat=ierror) rtemp
            if(field_label(i).gt.0)then
               unit_field(field_label(i),nunit) = rtemp
            endif
         enddo
         
         read(10, *, iostat=ierror) ctemp
      endif
   enddo
      
   close(10)
   
   call enter_exit(sub_name,2)
   
 end subroutine import_terminal
 
end module imports
