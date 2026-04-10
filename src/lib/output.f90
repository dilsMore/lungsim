subroutine exgauss(filename,groupname)
    character :: filename*(*),groupname*(*)
    !     Local Variables
    integer :: ifile=10,n,ne,ng,nj
    real(dp) :: xe(nsm,20),xi(3),ze(nsm,nhm),ze_nj(nsm),z(3)
    character :: readfile*150
    if(index(filename, ".exnode")>0)then ! full filename is given
       readfile = trim(filename)
    else! append correct extension
       readfile = trim(filename)//'.exnode'
    endif
    open(ifile, file=readfile, status='replace')
    !**   write the group name
    write(ifile,'( '' Group name: '',A)') trim(groupname)
    write(ifile,'('' #Fields=5'')')
    write(ifile,'('' 1) coordinates, coordinate,'',' &
         //''' rectangular cartesian, #Components=3'')')
    write(ifile,'(''   x. Value index=1, #Derivatives= 0'')')
    write(ifile,'(''   y. Value index=2, #Derivatives= 0'')')
    write(ifile,'(''   z. Value index=3, #Derivatives= 0'')')
    write(ifile,'(''  2) expansion_ratio, field,'',' &
         //''' rectangular cartesian, #Components=1'')')
    write(ifile,'(''   1. Value index=4, #Derivatives= 0'')')
    write(ifile,'(''  3) hydro_stress, field,'',' &
         //''' rectangular cartesian, #Components=1'')')
    write(ifile,'(''   1. Value index=5, #Derivatives= 0'')')
    write(ifile,'(''  4) compliance, field,'',' &
         //''' rectangular cartesian, #Components=1'')')
    write(ifile,'(''   1. Value index=6, #Derivatives= 0'')')
    write(ifile,'(''  5) material_label, field,'',' &
         //''' rectangular cartesian, #Components=1'')')
    write(ifile,'(''   1. Value index=7, #Derivatives= 0'')')
    n = 0
    do ne = 1,tissue_num_elems
       call xpxe(nb_lung,ne,xe)
       call zpze(nb_lung,ne,ze)
       do ng = 1,ngt(nb_lung)
          xi(1:3) = xig(1:3,ng,nb_lung)
          do nj = 1,3
             ze_nj(:) = ze(:,nj)
             z(nj)= PXI(nb_lung,1,xi,ze_nj)
          enddo
          n = n + 1
          write(ifile,'(1x,''Node: '',i12)') n
          write(ifile,'(2x, 3(f16.6))') z(1:3)
          write(ifile,'(2x, f16.6 )') yg(1,ng,ne)
          write(ifile,'(2x, f16.6 )') yg(2,ng,ne)/98.0665_dp
          write(ifile,'(2x, f16.6 )') yg(3,ng,ne)
          write(ifile,'(2x, f16.6 )') material_at_gp(5,ng,ne)
       enddo
    enddo                     !nolist (np)
    close(ifile)
  end subroutine exgauss

!#######################################################################

!subroutine export_1d_elem_field(ne_field, EXELEMFILE, group_name, field_name)
!    !!! TJ - output for airway unit
!    !!! Check list: Y/N
!
!    ! Element number (ne) - Y
!    ! xyz coordinates of the mid-point of the element (calculate from nodes) - Y
!    ! distance from entrance of model to the mid-point (calculate for each element) - Y (Not sure)
!    ! generation (elem_ordrs(1,ne)) - Y
!    ! Horsfield order (elem_ordrs(2,ne)) - Y
!    ! A label for the lobe that the element is in (calculate based on parent) - N
!    ! Length (elem_field(ne_length,ne)) - Y
!    ! radius (elem_field(ne_radius,ne)) - Y
!    ! flow rate (elem_field(ne_flow,ne)) - Y
!    ! mass in the lumen (elem_field(ne_mass,ne)) - Y
!    ! deposition by each mechanism - Y
!
!    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EXPORT_1D_ELEM_FIELD" :: EXPORT_1D_ELEM_FIELD
!
!!!! Parameters
!    integer, intent(in) :: ne_field
!    character(len=MAX_FILENAME_LEN), intent(in) :: EXELEMFILE
!    character(len=MAX_STRING_LEN), intent(in) :: field_name
!    character(len=MAX_STRING_LEN), intent(in) :: group_name
!
!!!! Local Variables
!    integer :: len_end,ne
!
!    integer :: np0
!    real(dp) :: midpoint
!
!    logical :: CHANGED
!    character(len=300) :: writefile
!    character(len=60) :: sub_name
!
!    ! --------------------------------------------------------------------------
!
!    sub_name = 'export_1d_elem_field'
!    call enter_exit(sub_name,1)
!
!    if(index(EXELEMFILE, ".exelem")> 0) then !full filename is given
!       writefile = EXELEMFILE
!    else ! need to append the correct filename extension
!       writefile = trim(EXELEMFILE)//'.exelem'
!    endif
!
!    open(10, file=writefile, status='replace')
!
!    len_end=len_trim(group_name)
!    !**     write the group name
!    write(10,'( '' Group name: '',A)') group_name(:len_end)
!    !**         write the elements
!    write(10,'( '' Shape.  Dimension=1'' )')
!    CHANGED=.TRUE. !initialise to force output of element information
!    len_end=len_trim(field_name)
!
!    ! set
!    np0 = 1 ! node number 1
!
!    do ne=1,num_elems
!       if(ne>1) THEN
!          CHANGED=.FALSE.
!       endif
!       if(CHANGED)THEN
!          write(10,'( '' #Scale factor sets=0'' )')
!          write(10,'( '' #Nodes= 0'' )')
!          write(10,'( '' #Fields= 1'' )')
!          write(10,'( '' 1)'',A,'', field, rectangular cartesian, #Components=1'')')&
!               field_name(:len_end)
!          write(10,'( ''  '',A,''.  l.Lagrange, no modify, grid based.'')') &
!               field_name(:len_end)
!          write(10,'( ''  #xi1=1'')')
!       endif
!
!       write(10,'(1X,''Element: '',I12,'' 0 0'' )') ne
!       !**               write the nodes
!       write(10,'(3X,''Nodes:'' )')
!       write(10,'(4X,2(1X,I12))') elem_nodes(1,ne),elem_nodes(2,ne)
!
!       !**               write the midpoints
!       write(10,'(3X,''Midpoint coordinate:'' )')
!       !write(10,'(4X,2(1X,I12))') (node_xyz(:, np1) + node_xyz(:, np2))/2
!       ! TJ - store midpoint coordinate
!       do nj=1,3
!            midpoint = (node_xyz(nj,np1) + node_xyz(nj, np2))/2
!            write(10,'(2X,4(1X,F12.6))') (node_xyz(nj,np1) + node_xyz(nj, np2))/2     !Coordinates
!       enddo !njj2
!
!       ! distance from entrance of model to the mid-point (calculate for each element)
!       !**               write the distance from entrance
!       write(10,'(3X,''Distance from Midpoint to entrance:'' )')
!       do nj=1,3
!           d = sqrt(sum(midpoint(nj, np)+node_xyz(nj, np0))**2)
!           write(10,'(2X,4(1X,F12.6))')  sqrt(sum(midpoint(nj, np) - node_xyz(nj, np0))**2) !distance
!       enddo !njj2
!
!
!       !**               write the generation
!       write(10,'(3X,''Generation:'' )')
!       write(10,'(4X,2(1X,E12.5))') elem_ordrs(1,ne)
!       !**               write the Horsfield order
!       write(10,'(3X,''Horsfield order:'' )')
!       write(10,'(4X,2(1X,E12.5))') elem_ordrs(2,ne)
!
!       ! A label for the lobe that the element is in (calculate based on parent)
!       ! Set default lung = 'Left'
!       ! set default lobe = 'superior'
!       !.flag. = true
!
!       !! suppose np = (x,y,z), where x-coordi defines the lung, z-coordi defines the lobe
!       !if (midpoint(1, np) .gt. 0) then direction = 'right' ! to find which lung
!
!       ! to find which lobe, set the range
!
!       if (midpoint(3, np) .lt. 0) then
!           lobe = 'right'
!       else
!       end if
!
!       !**               write the length
!       write(10,'(3X,''Length:'' )')
!       write(10,'(4X,2(1X,E12.5))') elem_ordrs(2,ne)
!
!       !**               write the radius
!       write(10,'(3X,''Radius:'' )')
!       write(10,'(4X,2(1X,E12.5))') elem_field(ne_radius,ne)   !
!       !**               write the flow rate
!       write(10,'(3X,''Flow rate:'' )')
!       write(10,'(4X,2(1X,E12.5))') elem_field(ne_flow,ne)   !
!       !**               write the mass in the lumen
!       write(10,'(3X,''Mass in lumen:'' )')
!       write(10,'(4X,2(1X,E12.5))') elem_field(ne_mass,ne)   !
!
!       !! deposition by each mechanism
!       ! bronchial dep mass by diffusion
!       write(10,'('' 4) nj_loss_dif, field, rectangular cartesian, #Components=1'')')
!       write(10,'(2X,4(1X,F12.6))') sum(node_field(nj_loss_dif,1:num_nodes))
!
!       ! bronchial dep mass by sedimentation
!       write(10,'('' 4) nj_loss_sed, field, rectangular cartesian, #Components=1'')')
!       write(10,'(2X,4(1X,F12.6))') sum(node_field(nj_loss_sed,1:num_nodes))
!
!       ! bronchial dep mass by impaction
!       write(10,'('' 4) nj_loss_imp, field, rectangular cartesian, #Components=1'')')
!       write(10,'(2X,4(1X,F12.6))')sum(node_field(nj_loss_imp,1:num_nodes))
!    enddo !no_nelist (ne)
!    close(10)
!
!  end subroutine export_1d_elem_field
!
!
!! ####################################################################################
!subroutine export_terminal_solution(EXNODEFILE, name)
!
!    !!! TJ - output for acinus unit
!    !!! Check list: Y/N
!    ! unit number - Y
!    ! xyz coordinates of the unit (this will be the end node location of the terminal element) - Y
!
!    ! distance from entrance - Y
!    ! A label for lobe the acinus is in - N
!
!    ! Initial volume of acinus - Y
!    ! Flow to the acinus - Y
!    ! Mass in the acinus lumen - Y
!    ! Initial volume of acinus - Y
!    ! Deposition in the acinus by each mechanism:
!    ! TJ - Merryn do we want the alveolar deposition for all acini or just the current single one?
!    ! alveolar dep mass by diffusion - Y.
!    ! alveolar dep mass by sedimentation - Y.
!
!    !!! Parameters
!    character(len=MAX_FILENAME_LEN),intent(in) :: EXNODEFILE
!    character(len=MAX_STRING_LEN),intent(in) :: name
!
!    !!! Local Variables
!    integer :: len_end,ne,nj,NOLIST,np,np_last,VALUE_INDEX
!    integer :: np0 ! initial node
!    logical :: FIRST_NODE
!
!    len_end=len_trim(name)
!    if(num_units.GT.0) THEN
!       open(10, file=EXNODEFILE, status='replace')
!       !**     write the group name
!       write(10,'( '' Group name: '',A)') name(:len_end)
!       VALUE_INDEX = 1
!       select case (model_type)
!        case ('particle_transport')
!          write(10,'( '' #Fields=4'' )')
!          write(10,'('' 1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
!          do nj=1,3
!             if(nj.eq.1) write(10,'(2X,''x.  '')',advance="no")
!             if(nj.eq.2) write(10,'(2X,''y.  '')',advance="no")
!             if(nj.eq.3) write(10,'(2X,''z.  '')',advance="no")
!             write(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
!             VALUE_INDEX=VALUE_INDEX+1
!          enddo
!
!          ! distance from entrance
!          write(10,'('' 2) Distance from entrance, field, rectangular cartesian, #Components=1'')')
!          write(10,'(2X,''1.  '')',advance="no")
!          write(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
!
!          ! A label for lobe the acinus is in
!!          write(10,'('' 2) Lobe, field, rectangular cartesian, #Components=1'')')
!!          write(10,'(2X,''1.  '')',advance="no")
!!          write(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
!
!          ! Initial volume of acinus
!          write(10,'('' 2) initial_volume, field, rectangular cartesian, #Components=1'')')
!          write(10,'(2X,''1.  '')',advance="no")
!          write(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
!
!          !Time averaged flow rate
!          write(10,'('' 2) average_flow, field, rectangular cartesian, #Components=1'')')
!          write(10,'(2X,''1.  '')',advance="no")
!          write(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
!          VALUE_INDEX=VALUE_INDEX+1
!
!          !unit mass
!          write(10,'('' 3) unit_mass, field, rectangular cartesian, #Components=1'')')
!          write(10,'(2X,''1.  '')',advance="no")
!          write(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
!          VALUE_INDEX=VALUE_INDEX+1
!
!          ! alveolar dep mass by diffusion
!          write(10,'('' 4) nu_loss_dif, field, rectangular cartesian, #Components=1'')')
!          write(10,'(2X,''1.  '')',advance="no")
!          write(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
!
!          ! alveolar dep mass by sedimentation
!          write(10,'('' 4) nu_loss_sed, field, rectangular cartesian, #Components=1'')')
!          write(10,'(2X,''1.  '')',advance="no")
!          write(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
!
!       end select
!
!    !*** Exporting Terminal Solution
!       do nolist = 1,num_units
!          ne = units(nolist)
!          np = elem_nodes(2,ne)
!          !***      write the node
!          write(10,'(1X,''Node: '',I12)') np
!          do nj=1,3
!             write(10,'(2X,4(1X,F12.6))') (node_xyz(nj,np))      !Coordinates
!          enddo !njj2
!          select case (model_type)
!            case ('particle_transport')
!             write(10,'(2X,4(1X,F12.6))') sqrt(sum(node_xyz(nj, np) - node_xyz(nj, np0))**2) ! distance from entrance
!             !write(10,'(2X,4(1X,F12.6))')   ! space left for Lobe
!             write(10,'(2X,4(1X,F12.6))') unit_field(nu_vol,nolist) ! Initial volume of acinus
!             write(10,'(2X,4(1X,F12.6))') (unit_field(nu_vdot0,nolist))   ! flow
!             write(10,'(2X,4(1X,F12.6))') (unit_field(nu_vol,nolist)*unit_field(nu_conc1,nolist)) ! unit_mass
!             ! TJ - we might want to change this if for the single acinus?
!             write(10,'(2X,4(1X,F12.6))') sum(unit_field(nu_loss_dif,1:num_units)) ! alveolar dep mass by dif
!             write(10,'(2X,4(1X,F12.6))') sum(unit_field(nu_loss_sed,1:num_units)) ! alveolar dep mass by sed
!          end select
!       enddo
!    endif
!end subroutine export_terminal_solution