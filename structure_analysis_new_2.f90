program structure_analysis
    implicit none
    
    integer, parameter ::  mxstep=5,nbins=200,endnbins=100,rgbin=250,natoms_peo=78,totresid=160,cpol1=1,cpol4=4,endatm_1=16,endatm_2=73
    integer::first_peo_atmid,rel_atmid
    real :: rbinsize,rgbinsize,endbinsize, x, y, z,hist_rdf_all(nbins),hist_rg(rgbin),hist_enddist(endnbins)
    real::x_dist,y_dist,z_dist,x_r,y_r,z_r,rr,total_mass,end_to_end_dist,rg_av,rdf_val,endtoend_av
    real :: dx,dy,dz,lx,ly,lz,rg_mol,r_lower,r_upper,r_center,V_shell,normfac,number_density
    real, allocatable ::  coord_wrap(:,:,:),coord_unwrap(:,:,:),ref_atmty(:),x_com(:),y_com(:),z_com(:),mass_peo(:),rg(:),endtoend(:)
    integer :: i,j,k,n_atoms,q,r,a,m,n,atm,step, atmty, atmid, ok,ix,iy,iz,resid,count_pair
    character(len=100) :: filename
    character(len=100) :: line
    real, parameter :: PI=4.D0*DATAN(1.D0)
    LOGICAL::pol_firstatm_found
    
    ! Initialize histogram bins
    rbinsize = 0.1 ! rdf bin width
    rgbinsize=0.1    
    endbinsize=0.5    ! 2 and 25
    count_pair=0
    allocate(mass_peo(7))
    mass_peo=[12.011,12.011,12.011,12.011,1.008,1.008,15.9994]
    
    filename="/home/vishnu/Documents/Vishnu/New_Build/polymer/peo/test/new_test/wrap_mod.lammpstrj"
    open(unit=10, file=filename, status='old')
    open(unit=11, file="/home/vishnu/Documents/Vishnu/New_Build/polymer/peo/test/new_test/rg_test.dat")
    open(unit=12, file="/home/vishnu/Documents/Vishnu/New_Build/polymer/peo/test/new_test/endtoend_test.dat")
    open(unit=13, file="/home/vishnu/Documents/Vishnu/New_Build/polymer/peo/test/new_test/totrdf_test.dat")
    step = 0 

    pol_firstatm_found= .false. 
    do
     ! Read the header to get the number of atoms
        if (step == 0) then
            do i = 1, 3
                read(10,*)line
                print*,line
            end do
            read(10,*) n_atoms
            print*,n_atoms,"natoms"
            ! Allocate arrays for storing coordinates
            allocate(coord_wrap(3, natoms_peo, totresid), coord_unwrap(3, natoms_peo, totresid),rg(totresid),ref_atmty(natoms_peo),endtoend(totresid))
            allocate(x_com(totresid),y_com(totresid),z_com(totresid))
            read(10,*)
            read(10,*)x,y
            lx=y-x 
            read(10,*)x,y
            ly=y-x 
            read(10,*)x,y
            lz=y-x 
            read(10,*)
            print*,lx,ly,lz 
            
        
        else 
            ! Skip header of the trajectory file
            read(10, *, iostat=ok) line
            if (ok < 0 .or. step .ge. mxstep) exit ! Exit loop if end of file is reached
            
            do i = 1, 8
                read(10,*)
            end do
            !print*,step,"step"
        end if

        ! print*,'polatom_no_outside',first_peo_atmid

        if (pol_firstatm_found==.false.) then
            first_peo_atmid = huge(0)
            do i=1,n_atoms
                read(10,*) atmty,atmid,x,y,z,ix,iy,iz,resid 
                if (atmty .le. 7) then 
                    if (atmid < first_peo_atmid) first_peo_atmid = atmid
                end if
                ! print*,"first peo atmid", first_peo_atmid
            end do
        end if
        rewind(10)
        ! Skip to atom section again

        do i = 1, 9
            read(10,*)
        end do

        do i = 1, n_atoms
            read(10,*) atmty,atmid,x,y,z,ix,iy,iz,resid    
            if (atmty .le. 7) then             ! atom type of PEO atoms 

                pol_firstatm_found= .True.
                ! print*,'inside if polatomfound',pol_firstatm_found
                ! print*,"atmid =",atmid
                
                rel_atmid= atmid-(first_peo_atmid-1)
                !  print*,'first_peo_atmid',first_peo_atmid
                r = MOD(rel_atmid - 1, natoms_peo) + 1
                q = (rel_atmid - 1) / natoms_peo + 1
                ! print*,'r and q',r,q
                if (q .ne. resid) then
                    print*,'q,resid',q, resid
                end if
                ref_atmty(r)=atmty
                coord_wrap(1, r, q) = x
                coord_wrap(2, r, q) = y
                coord_wrap(3, r, q) = z
                print*,'atmid,x',atmid,x

                coord_unwrap(1, r, q) = x+ix*lx
                coord_unwrap(2, r, q) = y+iy*ly
                coord_unwrap(3, r, q) = z+iz*lz
                print*,'unwrap x',atmid,coord_unwrap(1, r, q)
            end if

        end do
        step = step + 1
        do 42 i=1,totresid
            ! print*,'resid=',i
            x_com(i) = 0.0d0
            y_com(i) = 0.0d0
            z_com(i) = 0.0d0
            total_mass=0.d0
            rg_mol=0.0d0

            !!!!!!!!!!!!!!!!! END TO END DISTANCE !!!!!!!!!!!!!!!!!!!!!!!!!!!
            x_dist=coord_unwrap(1,endatm_1,i)-coord_unwrap(1,endatm_2,i)
            y_dist=coord_unwrap(2,endatm_1,i)-coord_unwrap(2,endatm_2,i)
            z_dist=coord_unwrap(3,endatm_1,i)-coord_unwrap(3,endatm_2,i)
            end_to_end_dist=sqrt(x_dist**2+y_dist**2+z_dist**2)
            endtoend(i)=endtoend(i)+end_to_end_dist
            m=1+floor(end_to_end_dist/endbinsize)
            hist_enddist(m)=hist_enddist(m)+1

            !!!!!!!!!!!!!!!!! END TO END DISTANCE !!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            do 43 j=1,natoms_peo    !looping through atoms in a single polymer
                ! print*,'aid=',j
                ! print*,'i',i
                atm=ref_atmty(j)
                x_com(i) = x_com(i)+(coord_unwrap(1,j,i)*mass_peo(atm))
                y_com(i) = y_com(i)+(coord_unwrap(2,j,i)*mass_peo(atm))
                z_com(i) = z_com(i)+(coord_unwrap(3,j,i)*mass_peo(atm))
                    
                total_mass = total_mass + mass_peo(atm)


                !!!!!!!!!!!!!!!!! RDF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                do k=j+1,natoms_peo
                    x_r=coord_wrap(1,j,i)-coord_wrap(1,k,i)
                    y_r=coord_wrap(2,j,i)-coord_wrap(2,k,i)
                    z_r=coord_wrap(3,j,i)-coord_wrap(3,k,i)
                    rr= sqrt(x_r**2+y_r**2+z_r**2)
                    !print*,'rr',rr,i,j,k,coord_unwrap(1,j,i),coord_unwrap(1,k,i)
                    a=1+floor(rr/rbinsize)
                    if(rr<15) then 
                    !print*,"1"
                    count_pair=count_pair+1
                    hist_rdf_all(a)= hist_rdf_all(a) + 1
                    end if
                end do
                !!!!!!!!!!!!!!!!! RDF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            43 end do
            x_com(i) = x_com(i) / total_mass
            y_com(i) = y_com(i) / total_mass
            z_com(i) = z_com(i) / total_mass
            ! print*,'x_com(i)',x_com(i)
            ! print*,'total_mass',total_mass



            !!!!!!!!!!!!!!!!! RG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do j = 1,natoms_peo
                dx = coord_unwrap(1,j,i) - x_com(i)
                dy = coord_unwrap(2,j,i) - y_com(i)
                dz = coord_unwrap(3,j,i) - z_com(i)
                rg_mol = rg_mol + (mass_peo(ref_atmty(j))*(dx*dx + dy*dy + dz*dz))/total_mass
                ! print*,'rg_mol',i,j,rg_mol
            end do
            rg_mol=sqrt(rg_mol)
            ! print*,'rg(i) bi',rg(i)
            rg(i) = rg(i)+rg_mol
            n=1+floor(rg_mol/rgbinsize)
            hist_rg(n)=hist_rg(n)+1
            ! print*,rg(i)
            !!!!!!!!!!!!!!!!! RG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        42 end do
    end do

    print*,'step',step
    print*,'natoms',n_atoms
    print*,'count_pair',count_pair

    do j =1,totresid       !!!!!!!!!! Average Rg over frames !!!!!!!!!!!!!!!!!!
        !print*,"radius of gyration for each molecule",rg(j)/step
        rg_av=rg_av+rg(j)/step
        endtoend_av=endtoend_av+(endtoend(j)/step)
        !print*,"end to end distance for each molecule",endtoend(j)/step
    end do
    print*,'Average radius of gyration=',rg_av/totresid
    print*,'Average endtoend distance=',endtoend_av/totresid

     number_density = dble(totresid*78) / (lx * ly * lz)

    do n = 1, rgbin
        r_lower = (n - 1) * rgbinsize
        r_upper = n * rgbinsize
        r_center = 0.5d0 * (r_lower + r_upper)      
        
        hist_rg(n)= hist_rg(n) / dble(totresid * step)
        write(11, '(2F12.4)') r_center,hist_rg(n)
     end do
    close(11)

    do n = 1, endnbins
        r_lower = (n - 1) * endbinsize
        r_upper = n * endbinsize
        r_center = 0.5d0 * (r_lower + r_upper)      
        
        hist_enddist(n)=hist_enddist(n)/dble(totresid * step)
        write(12, '(15F12.4)') r_center,hist_enddist(n)
    end do
    close(12)
     
    do n = 1, nbins
        r_lower = (n - 1) * rbinsize
        r_upper = n * rbinsize
        r_center = 0.5d0 * (r_lower + r_upper)
        V_shell = (4.0d0 / 3.0d0) * PI * (r_upper**3 - r_lower**3)
        !normfac= number_density * dble(totresid * step) * V_shell
        normfac= dble(totresid * step) * V_shell
        rdf_val = hist_rdf_all(n) / normfac
        !normfac = dble(count_pair) * V_shell / (boxx * boxy * boxz)
        !normfac=(4.0/3.0)*PI*(3*(n**2)-3*n+1)*(binsize**3)*count_pair/(2*boxx*boxy*boxz)
        write(13, '(15F12.4)') r_center,rdf_val
    end do
    close(13)


    
    close(10)
end program structure_analysis

