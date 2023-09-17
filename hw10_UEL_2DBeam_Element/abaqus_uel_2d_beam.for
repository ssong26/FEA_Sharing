!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 1D Timosike beam
!
!    The file also contains the following subrouines:
!          abq_UEL_2D_integrationpoints           - defines integration points for 2D continuum elements
!          abq_UEL_2D_shapefunctions              - defines shape functions for 2D continuum elements
!          abq_UEL_1D_integrationpoints(n_points, n_nodes, xi, w)  = defines integration points for 1D line integral
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3     LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
    !
      INCLUDE 'ABA_PARAM.INC'
    !
    !
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1   SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2   DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3   JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4   PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

    !
    !       Variables that must be computed in this routine
    !       RHS(i)                     Right hand side vector.  In EN234_FEA the dimensions are always RHS(MLVARX,1)
    !       AMATRX(i,j)                Stiffness matrix d RHS(i)/ d DU(j)
    !       SVARS(1:NSVARS)            Element state variables.  Must be updated in this routine
    !       ENERGY(1:8)
    !                                  Energy(1) Kinetic Energy
    !                                  Energy(2) Elastic Strain Energy
    !                                  Energy(3) Creep Dissipation
    !                                  Energy(4) Plastic Dissipation
    !                                  Energy(5) Viscous Dissipation
    !                                  Energy(6) Artificial strain energy
    !                                  Energy(7) Electrostatic energy
    !                                  Energy(8) Incremental work done by loads applied to the element
    !       PNEWDT                     Allows user to control ABAQUS time increments.
    !                                  If PNEWDT<1 then time step is abandoned and computation is restarted with
    !                                  a time increment equal to PNEWDT*DTIME
    !                                  If PNEWDT>1 ABAQUS may increase the time increment by a factor PNEWDT
    !
    !       Variables provided for information
    !       NDOFEL                     Total # DOF for the element
    !       NRHS                       Dimension variable
    !       NSVARS                     Total # element state variables
    !       PROPS(1:NPROPS)            User-specified properties of the element
    !       NPROPS                     No. properties
    !       JPROPS(1:NJPROPS)          Integer valued user specified properties for the element
    !       NJPROPS                    No. integer valued properties
    !       COORDS(i,N)                ith coordinate of Nth node on element
    !       MCRD                       Maximum of (# coords,minimum of (3,#DOF)) on any node
    !       U                          Vector of DOF at the end of the increment
    !       DU                         Vector of DOF increments
    !       V                          Vector of velocities (defined only for implicit dynamics)
    !       A                          Vector of accelerations (defined only for implicit dynamics)
    !       JTYPE                      Integer identifying element type (the number n in the Un specification in the input file)
    !       TIME(1:2)                  TIME(1)   Current value of step time
    !                                  TIME(2)   Total time
    !       DTIME                      Time increment
    !       KSTEP                      Current step number (always 1 in EN234_FEA)
    !       KINC                       Increment number
    !       JELEM                      User assigned element number in ABAQUS (internally assigned in EN234_FEA)
    !       PARAMS(1:3)                Time increment parameters alpha, beta, gamma for implicit dynamics
    !       NDLOAD                     Number of user-defined distributed loads defined for this element
    !       JDLTYP(1:NDLOAD)           Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
    !       ADLMAG(1:NDLOAD)           Distributed load magnitudes
    !       DDLMAG(1:NDLOAD)           Increment in distributed load magnitudes
    !       PREDEF(1:2,1:NPREDF,1:NNODE)   Predefined fields.
    !       PREDEF(1,...)              Value of predefined field
    !       PREDEF(2,...)              Increment in predefined field
    !       PREDEF(1:2,1,k)            Value of temperature/temperature increment at kth node
    !       PREDEF(1:2,2:NPREDF,k)     Value of user defined field/field increment at kth node (not used in EN234FEA)
    !       NPREDF                     Number of predefined fields (1 for en234FEA)
    !       LFLAGS                     Control variable
    !       LFLAGS(1)                  Defines procedure type
    !       LFLAGS(2)                  0 => small displacement analysis  1 => Large displacement (NLGEOM option)
    !       LFLAGS(3)                   1 => Subroutine must return both RHS and AMATRX (always true in EN234FEA)
    !                                   2 => Subroutine must return stiffness AMATRX = -dF/du
    !                                   3 => Subroutine must return daming matrix AMATRX = -dF/dudot
    !                                   4 => Subroutine must return mass matrix AMATRX = -dF/duddot
    !                                   5 => Define the RHS only
    !                                   6 => Define the mass matrix for the initial acceleration calculation
    !                                   100 => Define perturbation quantities for output
    !       LFLAGS(4)                   0 => General step   1 => linear perturbation step
    !       LFLAGS(5)                   0 => current approximation to solution based on Newton correction; 1 => based on extrapolation
    !       MLVARX                      Dimension variable (equal to NDOFEL in EN234FEA)
    !       PERIOD                      Time period of the current step
    !
    !
    ! Local Variables
      integer      :: i,j,n_points,kint, nfacenodes
      integer      :: face_node_list(3)                       ! List of nodes on an element face
    !
      double precision  ::  xi(2,9)                          ! Area integration points
      double precision  ::  w(5)                             ! Area integration weights
      double precision  ::  N(9)                             ! 2D shape functions
      double precision  ::  dNdxi(9,2)                       ! 2D shape function derivatives
      double precision  ::  dNdx(9,2)                        ! Spatial derivatives
      double precision  ::  dxdxi(2,2)                       ! Derivative of spatial coords wrt normalized coords

    !   Variables below are for computing integrals over element faces
      double precision  ::  face_coords(2,3)                  ! Coords of nodes on an element face
      double precision  ::  xi1(5)                            ! 1D integration points
      double precision  ::  w1(5)                             ! Integration weights
      double precision  ::  N1(3)                             ! 1D shape functions
      double precision  ::  dN1dxi(3)                         ! 1D shape function derivatives
      double precision  ::  norm(2)                           ! Normal to an element face
      double precision  ::  dxdxi1(2)                         ! Derivative of 1D spatial coord wrt normalized areal coord
      double precision  ::  xi2(2,5)                          ! Integration points in 2D.  
      double precision  ::  dxdxi2(2)
    !
      double precision  ::  strain_star(2)                    ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
      double precision  ::  stress_star(2)                    ! Stress vector contains [s11, s22, s33, s12, s13, s23]
      double precision  ::  D(2,2)                            ! stress = D*(strain)
      double precision  ::  B(3,8)                            ! strain = B*(dof_total)
      double precision  ::  T(8,6)                            ! T matrix
      double precision  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
      double precision  ::  E, xnu                            ! Material properties
      double precision  ::  G                                 ! Shear modulus
      double precision  ::  thickness, width                  ! Beam Properties
      double precision  ::  utemp(22)                         ! Temporary displacement (for incompatible mode elements)
      double precision  ::  cos_1,sin_1                       ! Cos_1 = cos_theta1, Sin_1 = sin_theta1
      double precision  ::  cos_2,sin_2                       ! Cos_2 = cos_theta2, Sin_2 = sin_theta2
      double precision  ::  Length                            ! Store the local length
      double precision  ::  coord_slave(2,4)                  ! coord_slave is to store the coordinate for the slave points
      double precision  ::  estar(2,2)                        ! Laminar basis vectors
      double precision  ::  R(2,3)                            ! Rotation tensor
      double precision  ::  RBT(2,6)                          ! RBT = matmul(R,B,T)
      double precision  ::  axial_force                       ! axial force
      double precision  ::  shear_force                       ! shear force
      double precision  ::  bending_moment                    ! bending moment

    !
    !     Example ABAQUS UEL implementing 2D linear elastic elements
    !     Includes option for incompatible mode elements
    !     El props are:
    
    !     PROPS(1)         thickness
    !     PROPS(2)         width
    !     PROPS(3)         Young's modulus
    !     PROPS(4)         Poisson's ratio
    !
      ! ======================================================================================
      ! =============================Store the Material Properties============================
      !
      thickness = PROPS(1)
      width = PROPS(2)
      E = PROPS(3)
      xnu = PROPS(4)
      
      ! ======================================================================================
      ! =============================Store the Structural Properties==========================
      !
      ! The length of the local coordinate
      Length = sqrt((coords(1,2)-coords(1,1))**2.d0 + ! coords(x1:x3,1:NNODE)
     1(coords(2,2)-coords(2,1))**2.d0)
      
      ! The angular properties.
      cos_1 = (coords(1,2)-coords(1,1))/Length
      sin_1 = (coords(2,2)-coords(2,1))/Length
      cos_2 = (coords(1,2)-coords(1,1))/Length
      sin_2 = (coords(2,2)-coords(2,1))/Length
      
      ! To store the slave coordinate, for calculation in the intergration loop.
      ! Point 1
      coord_slave(1,1) = coords(1,1) + thickness/2.d0*sin_1
      coord_slave(2,1) = coords(2,1) - thickness/2.d0*cos_1
      ! Point 2
      coord_slave(1,2) = coords(1,2) + thickness/2.d0*sin_2
      coord_slave(2,2) = coords(2,2) - thickness/2.d0*cos_2
      ! Point 3
      coord_slave(1,3) = coords(1,2) - thickness/2.d0*sin_2
      coord_slave(2,3) = coords(2,2) + thickness/2.d0*cos_2
      ! Point 4
      coord_slave(1,4) = coords(1,1) - thickness/2.d0*sin_1
      coord_slave(2,4) = coords(2,1) + thickness/2.d0*cos_1
      
      ! In order to use a 5 points trapozoidal scheme
      n_points = 5

      ! ======================================================================================
      ! =============================Intergration over all points=============================
      
      ! Intergration points, we use 5 points trapozoidal scheme 
      call abq_UEL_1D_integrationpoints(n_points, NNODE, xi1, w) ! xi1 is the local points, w is the weight
      
      ! To use the five points in a 2-D condition
      xi2 = 0.d0
      xi2(2,1:5) = xi1(1:5)
      
      ! Initialize output parameters
      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
      
      ! Initialize state parameters
      axial_force = 0.d0
      shear_force = 0.d0
      bending_moment = 0.d0
      
      ! Initialize D matrix
      G = E/2.d0/(1+xnu)
      D = 0.d0
      D(1,1) = E
      D(2,2) = G
      
      ! Initialize T matrix
      T = 0.d0
      
      T(1,1) = 1.d0
      T(1,3) = coords(2,1) - coord_slave(2,1)
      
      T(2,2) = 1.d0
      T(2,3) = - (coords(1,1) - coord_slave(1,1))
      
      T(3,4) = 1.d0
      T(3,6) = coords(2,2) - coord_slave(2,2)
      
      T(4,5) = 1.d0
      T(4,6) = coords(1,2) - coord_slave(1,2)
      
      T(5,4) = 1.d0
      T(5,6) = coords(2,2) - coord_slave(2,3)
      
      T(6,5) = 1.d0
      T(6,6) = coords(1,2) - coord_slave(1,3)
      
      T(7,1) = 1.d0
      T(7,3) = coords(2,1) - coord_slave(2,4)
      
      T(8,2) = 1.d0
      T(8,3) = - (coords(1,1) - coord_slave(1,4))
      
      ! Initialize energy
      ENERGY(1:8) = 0.d0

      if (MLVARX<2*NNODE) then
        write(6,*) ' Error in abaqus UEL '
        write(6,*) ' Variable MLVARX must exceed 2*NNODE'
        write(6,*) ' MLVARX = ',MLVARX,' NNODE = ',NNODE
        stop
      endif
      
      ! Loop over all points
        do kint = 1, n_points
            
            ! Use the four nodes elements to determine B matrix
            call abq_UEL_2D_shapefunctions(xi2(1:2,kint),4,N,dNdxi)
            dxdxi = matmul(coord_slave(1:2,1:4),dNdxi(1:4,1:2))
            
            call abq_inverse_LU(dxdxi,dxidx,2)
            dNdx(1:4,1:2) = matmul(dNdxi(1:4,1:2),dxidx)
            
            ! Determine B matrix
            B = 0.d0
            B(1,1:2*4-1:2) = dNdx(1:4,1)
            B(2,2:2*4:2)   = dNdx(1:4,2)
            B(3,1:2*4-1:2) = dNdx(1:4,2)
            B(3,2:2*4:2)   = dNdx(1:4,1)
            
            ! determine estar
            estar(1:2,1) = dxdxi(1:2,1)/norm2(dxdxi(1:2,1))
            estar(1:2,2) = dxdxi(1:2,2)/norm2(dxdxi(1:2,2))
            
            ! To decide R matrix
            R = 0.d0
            
            R(1,1) = estar(1,1)**2.d0
            R(1,2) = estar(2,1)**2.d0
            R(1,3) = estar(1,1)*estar(2,1)
            
            R(2,1) = 2*estar(1,1)*estar(1,2)
            R(2,2) = 2*estar(2,1)*estar(2,2)
            R(2,3) = estar(1,1)*estar(2,2) + estar(1,2)*estar(2,1)
            
            RBT = matmul(matmul(R,B),T)
            ! To determine strain and stress
            strain_star = matmul(RBT,U(1:6))
            stress_star = matmul(D,strain_star)
            
            RHS(1:6,1) = RHS(1:6,1)
     1      - matmul(transpose(RBT),stress_star)*w(kint)*width*Length
     2      * thickness/2.d0

            AMATRX(1:6,1:6) = AMATRX(1:6,1:6)
     1      + matmul(matmul(transpose(RBT),D),RBT)*
     2      w(kint)*width*Length*thickness/2.d0

            ! Store state variables
            axial_force = axial_force
     1      + stress_star(1)*w(kint)*width*thickness/2.d0
            
            shear_force = shear_force
     1      + stress_star(2)*w(kint)*width*thickness/2.d0
            
            bending_moment = bending_moment
     1      + xi2(2,kint)*stress_star(1)*w(kint)*width*thickness**2.d0
     2      /4.d0
        end do
        
      ! To store the state viarable
      SVARS(1) = axial_force
      SVARS(2) = shear_force
      SVARS(3) = bending_moment
      
      ! Unknown
      PNEWDT = 1.d0          ! This leaves the timestep unchanged (ABAQUS will use its own algorithm to determine DTIME)
      
      return
      
      END SUBROUTINE UEL

      subroutine abq_UEL_2D_integrationpoints(n_points, n_nodes, xi, w)
    
      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes
    
      double precision, intent(out) :: xi(2,*)
      double precision, intent(out) :: w(*)
    
      integer :: i,j,k,n
    
      double precision :: cn,w1,w2,w11,w12,w22
    
    !         Defines integration points and weights for 2D continuum elements
    
      if ( n_points==1 ) then
        if ( n_nodes==4 .or. n_nodes==9 ) then    !     ---   4 or 9 noded quad
            xi(1, 1) = 0.D0
            xi(2, 1) = 0.D0
            w(1) = 4.D0
        else if ( n_nodes==3 .or. n_nodes==6 ) then !     ---   3 or 6 noded triangle
            xi(1, 1) = 1.D0/3.D0
            xi(2, 1) = 1.D0/3.D0
            w(1) = 1.D0/2.D0
        end if
      else if ( n_points==3 ) then
        xi(1, 1) = 0.5D0
        xi(2, 1) = 0.5D0
        w(1) = 1.D0/6.D0
        xi(1, 2) = 0.D0
        xi(2, 2) = 0.5D0
        w(2) = w(1)
        xi(1, 3) = 0.5D0
        xi(2, 3) = 0.D0
        w(3) = w(1)
      else if ( n_points==4 ) then
        if ( n_nodes==4 .or. n_nodes==8 .or. n_nodes==9 ) then
            !     2X2 GAUSS INTEGRATION POINTS FOR QUADRILATERAL
            !     43
            !     12
            cn = 0.5773502691896260D0
            xi(1, 1) = -cn
            xi(1, 2) = cn
            xi(1, 3) = cn
            xi(1, 4) = -cn
            xi(2, 1) = -cn
            xi(2, 2) = -cn
            xi(2, 3) = cn
            xi(2, 4) = cn
            w(1) = 1.D0
            w(2) = 1.D0
            w(3) = 1.D0
            w(4) = 1.D0
        else if ( n_nodes==3 .or. n_nodes==6 ) then
            !     xi integration points for triangle
            xi(1, 1) = 1.D0/3.D0
            xi(2, 1) = xi(1, 1)
            w(1) = -27.D0/96.D0
            xi(1, 2) = 0.6D0
            xi(2, 2) = 0.2D0
            w(2) = 25.D0/96.D0
            xi(1, 3) = 0.2D0
            xi(2, 3) = 0.6D0
            w(3) = w(2)
            xi(1, 4) = 0.2D0
            xi(2, 4) = 0.2D0
            w(4) = w(2)
        end if
    
      else if ( n_points==7 ) then
        ! Quintic integration for triangle
        xi(1,1) = 1.d0/3.d0
        xi(2,1) = xi(1,1)
        w(1) = 0.1125d0
        xi(1,2) = 0.0597158717d0
        xi(2,2) = 0.4701420641d0
        w(2) = 0.0661970763d0
        xi(1,3) = xi(2,2)
        xi(2,3) = xi(1,2)
        w(3) = w(2)
        xi(1,4) = xi(2,2)
        xi(2,4) = xi(2,2)
        w(4) = w(2)
        xi(1,5) = 0.7974269853d0
        xi(2,5) = 0.1012865073d0
        w(5) = 0.0629695902d0
        xi(1,6) = xi(2,5)
        xi(2,6) = xi(1,5)
        w(6) = w(5)
        xi(1,7) = xi(2,5)
        xi(2,7) = xi(2,5)
        w(7) = w(5)
      else if ( n_points==9 ) then
        !     3X3 GAUSS INTEGRATION POINTS
        !     789
        !     456
        !     123
        cn = 0.7745966692414830D0
        xi(1, 1) = -cn
        xi(1, 2) = 0.D0
        xi(1, 3) = cn
        xi(1, 4) = -cn
        xi(1, 5) = 0.D0
        xi(1, 6) = cn
        xi(1, 7) = -cn
        xi(1, 8) = 0.D0
        xi(1, 9) = cn
        xi(2, 1) = -cn
        xi(2, 2) = -cn
        xi(2, 3) = -cn
        xi(2, 4) = 0.D0
        xi(2, 5) = 0.D0
        xi(2, 6) = 0.D0
        xi(2, 7) = cn
        xi(2, 8) = cn
        xi(2, 9) = cn
        w1 = 0.5555555555555560D0
        w2 = 0.8888888888888890D0
        w11 = w1*w1
        w12 = w1*w2
        w22 = w2*w2
        w(1) = w11
        w(2) = w12
        w(3) = w11
        w(4) = w12
        w(5) = w22
        w(6) = w12
        w(7) = w11
        w(8) = w12
        w(9) = w11
      end if
    
      return
    
      end subroutine abq_UEL_2D_integrationpoints   
    
      subroutine abq_UEL_2D_shapefunctions(xi,n_nodes,f,df)
    
      implicit none
      integer, intent(in) :: n_nodes
    
      double precision, intent(in) :: xi(2)
      double precision, intent(out) :: f(*)
      double precision, intent(out) :: df(9,2)
      double precision g1, g2, g3, dg1, dg2, dg3
      double precision h1, h2, h3, dh1, dh2, dh3
      double precision z,dzdp, dzdq
    
            if ( n_nodes==3 ) then        !     SHAPE FUNCTIONS FOR 3 NODED TRIANGLE
                f(1) = xi(1)
                f(2) = xi(2)
                f(3) = 1.D0 - xi(1) - xi(2)
                df(1, 1) = 1.D0
                df(1, 2) = 0.D0
                df(2, 1) = 0.D0
                df(2, 2) = 1.D0
                df(3, 1) = -1.D0
                df(3, 2) = -1.D0
            else if ( n_nodes==4 ) then
                !     SHAPE FUNCTIONS FOR 4 NODED QUADRILATERAL
                !     43
                !     12
                g1 = 0.5D0*(1.D0 - xi(1))
                g2 = 0.5D0*(1.D0 + xi(1))
                h1 = 0.5D0*(1.D0 - xi(2))
                h2 = 0.5D0*(1.D0 + xi(2))
                f(1) = g1*h1
                f(2) = g2*h1
                f(3) = g2*h2
                f(4) = g1*h2
                dg1 = -0.5D0
                dg2 = 0.5D0
                dh1 = -0.5D0
                dh2 = 0.5D0
                df(1, 1) = dg1*h1
                df(2, 1) = dg2*h1
                df(3, 1) = dg2*h2
                df(4, 1) = dg1*h2
                df(1, 2) = g1*dh1
                df(2, 2) = g2*dh1
                df(3, 2) = g2*dh2
                df(4, 2) = g1*dh2
    
            else if ( n_nodes==6 ) then
    
                !     SHAPE FUNCTIONS FOR 6 NODED TRIANGLE
                !          3
    
                !       6      5
    
                !     1    4     2
    
                !     P = L1
                !     Q = L2
                !     Z = 1 - P - Q = L3
    
                z = 1.D0 - xi(1) - xi(2)
                f(1) = (2.D0*xi(1) - 1.D0)*xi(1)
                f(2) = (2.D0*xi(2) - 1.D0)*xi(2)
                f(3) = (2.D0*z - 1.D0)*z
                f(4) = 4.D0*xi(1)*xi(2)
                f(5) = 4.D0*xi(2)*z
                f(6) = 4.D0*xi(1)*z
                dzdp = -1.D0
                dzdq = -1.D0
                df(1, 1) = 4.D0*xi(1) - 1.D0
                df(2, 1) = 0.D0
                df(3, 1) = 4.D0*z*dzdp - dzdp
                df(4, 1) = 4.D0*xi(2)
                df(5, 1) = 4.D0*xi(2)*dzdp
                df(6, 1) = 4.D0*z + 4.D0*xi(1)*dzdp
                df(1, 2) = 0.D0
                df(2, 2) = 4.D0*xi(2) - 1.D0
                df(3, 2) = 4.D0*z*dzdq - dzdq
                df(4, 2) = 4.D0*xi(1)
                df(5, 2) = 4.D0*z + 4.D0*xi(2)*dzdq
                df(6, 2) = 4.D0*xi(1)*dzdq
    
            else if ( n_nodes==8 ) then
                !     SHAPE FUNCTIONS FOR 8 NODED SERENDIPITY ELEMENT
                 f(1) = -0.25*(1.-xi(1))*(1.-xi(2))*(1.+xi(1)+xi(2));
                 f(2) = 0.25*(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-1.);
                 f(3) = 0.25*(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-1.);
                 f(4) = 0.25*(1.-xi(1))*(1.+xi(2))*(xi(2)-xi(1)-1.);
                 f(5) = 0.5*(1.-xi(1)*xi(1))*(1.-xi(2));
                 f(6) = 0.5*(1.+xi(1))*(1.-xi(2)*xi(2));
                 f(7) = 0.5*(1.-xi(1)*xi(1))*(1.+xi(2));
                 f(8) = 0.5*(1.-xi(1))*(1.-xi(2)*xi(2));
                 df(1,1) = 0.25*(1.-xi(2))*(2.*xi(1)+xi(2));
                 df(1,2) = 0.25*(1.-xi(1))*(xi(1)+2.*xi(2));
                 df(2,1) = 0.25*(1.-xi(2))*(2.*xi(1)-xi(2));
                 df(2,2) = 0.25*(1.+xi(1))*(2.*xi(2)-xi(1));
                 df(3,1) = 0.25*(1.+xi(2))*(2.*xi(1)+xi(2));
                 df(3,2) = 0.25*(1.+xi(1))*(2.*xi(2)+xi(1));
                 df(4,1) = 0.25*(1.+xi(2))*(2.*xi(1)-xi(2));
                 df(4,2) = 0.25*(1.-xi(1))*(2.*xi(2)-xi(1));
                 df(5,1) = -xi(1)*(1.-xi(2));
                 df(5,2) = -0.5*(1.-xi(1)*xi(1));
                 df(6,1) = 0.5*(1.-xi(2)*xi(2));
                 df(6,2) = -(1.+xi(1))*xi(2);
                 df(7,1) = -xi(1)*(1.+xi(2));
                 df(7,2) = 0.5*(1.-xi(1)*xi(1));
                 df(8,1) = -0.5*(1.-xi(2)*xi(2));
                 df(8,2) = -(1.-xi(1))*xi(2);
            else if ( n_nodes==9 ) then
                !     SHAPE FUNCTIONS FOR 9 NODED LAGRANGIAN ELEMENT
                !     789
                !     456
                !     123
                g1 = -.5D0*xi(1)*(1.D0 - xi(1))
                g2 = (1.D0 - xi(1))*(1.D0 + xi(1))
                g3 = .5D0*xi(1)*(1.D0 + xi(1))
                h1 = -.5D0*xi(2)*(1.D0 - xi(2))
                h2 = (1.D0 - xi(2))*(1.D0 + xi(2))
                h3 = .5D0*xi(2)*(1.D0 + xi(2))
                dg1 = xi(1) - 0.5d0
                dg2 = -2.d0*xi(1)
                dg3 = xi(1) + 0.5d0
                dh1 = xi(2)-0.5d0
                dh2 = -2.d0*xi(2)
                dh3 = xi(2) + 0.5d0
                f(1) = g1*h1
                f(2) = g2*h1
                f(3) = g3*h1
                f(4) = g1*h2
                f(5) = g2*h2
                f(6) = g3*h2
                f(7) = g1*h3
                f(8) = g2*h3
                f(9) = g3*h3
                df(1,1) = dg1*h1
                df(1,2) = g1*dh1
                df(2,1) = dg2*h1
                df(2,2) = g2*dh1
                df(3,1) = dg3*h1
                df(3,2) = g3*dh1
                df(4,1) = dg1*h2
                df(4,2) = g1*dh2
                df(5,1) = dg2*h2
                df(5,2) = g2*dh2
                df(6,1) = dg3*h2
                df(6,2) = g3*dh2
                df(7,1) = dg1*h3
                df(7,2) = g1*dh3
                df(8,1) = dg2*h3
                df(8,2) = g2*dh3
                df(9,1) = dg3*h3
                df(9,2) = g3*dh3
            end if
    
      end subroutine abq_UEL_2D_shapefunctions
     
      subroutine abq_UEL_1D_integrationpoints(n_points, n_nodes, xi, w)
    
    
      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes
    
      double precision, intent(out) :: xi(*)
      double precision, intent(out) :: w(*)
    
      integer :: i,j,k,n
    
      double precision x1D(4), w1D(4)
    
    
    
      select case ( n_points )
        case (2)
            xi(1) = .5773502691896257D+00
            xi(2) = -.5773502691896257D+00
            w(1) = .1000000000000000D+01
            w(2) = .1000000000000000D+01
            return
        case (3)
            xi(1) = 0.7745966692414834D+00
            xi(2) = .0000000000000000D+00
            xi(3) = -.7745966692414834D+00
            w(1) = .5555555555555556D+00
            w(2) = .8888888888888888D+00
            w(3) = .5555555555555556D+00
            return
        case (4)
            xi(1) = .8611363115940526D+00
            xi(2) = .3399810435848563D+00
            xi(3) = -.3399810435848563D+00
            xi(4) = -.8611363115940526D+00
            w(1) = .3478548451374538D+00
            w(2) = .6521451548625461D+00
            w(3) = .6521451548625461D+00
            w(4) = .3478548451374538D+00
            return
        case (5)
            xi(1) = .9061798459386640D+00
            xi(2) = .5384693101056831D+00
            xi(3) = .0000000000000000D+00
            xi(4) = -.5384693101056831D+00
            xi(5) = -.9061798459386640D+00
            w(1) = .2369268850561891D+00
            w(2) = .4786286704993665D+00
            w(3) = .5688888888888889D+00
            w(4) = .4786286704993665D+00
            w(5) = .2369268850561891D+00
            return
        case (6)
            xi(1) = .9324695142031521D+00
            xi(2) = .6612093864662645D+00
            xi(3) = .2386191860831969D+00
            xi(4) = -.2386191860831969D+00
            xi(5) = -.6612093864662645D+00
            xi(6) = -.9324695142031521D+00
            w(1) = .1713244923791703D+00
            w(2) = .3607615730481386D+00
            w(3) = .4679139345726910D+00
            w(4) = .4679139345726910D+00
            w(5) = .3607615730481386D+00
            w(6) = .1713244923791703D+00
            return
        case DEFAULT
            write(6,*)'Error in subroutine abq_UEL_1D_integrationpoints'
            write(6,*) ' Invalid number of integration points for 1D'
            write(6,*) ' n_points must be between 1 and 6'
            stop
      end select
    
    
    
    
    
    
    
      end subroutine ABQ_UEL_1D_integrationpoints  
    
      subroutine abq_facenodes_2D(nelnodes,face,list,nfacenodes)
    
      implicit none
    
      integer, intent (in)      :: nelnodes
      integer, intent (in)      :: face
      integer, intent (out)     :: list(*)
      integer, intent (out)     :: nfacenodes
    !
    !        Subroutine to return list of nodes on an element face for standard 2D solid elements
    !
      integer :: i3(3)
      integer :: i4(4)
    
      i3(1:3) = [2,3,1]
      i4(1:4) = [2,3,4,1]
    
      if (nelnodes == 3) then
        nfacenodes = 2
        list(1) = face
        list(2) = i3(face)
      else if (nelnodes == 4) then
        nfacenodes = 2
        list(1) = face
        list(2) = i4(face)
      else if (nelnodes == 6) then
        nfacenodes = 3
        list(1) = face
        list(2) = i3(face)
        list(3) = face+3
      else if (nelnodes == 8) then
        nfacenodes = 3
        list(1) = face
        list(2) = i4(face)
        list(3) = face+4
      else if (nelnodes == 9) then
        nfacenodes = 3
        if (face==1) list(1:3) = (/1,3,2/)
        if (face==2) list(1:3) = (/3,9,6/)
        if (face==3) list(1:3) = (/9,7,8/)
        if (face==4) list(1:3) = (/7,1,4/)
      endif
    
      end subroutine abq_facenodes_2d
    
      subroutine abq_inverse_LU(Ain,A_inverse,n)  ! Compute the inverse of an arbitrary matrix by LU decomposition
    
        implicit none
    
        integer, intent(in)  :: n
    
        double precision, intent(in)    :: Ain(n,n)
        double precision, intent(out)   :: A_inverse(n,n)
    
        double precision :: A(n,n), L(n,n), U(n,n), b(n), d(n), x(n)
        double precision :: coeff
        integer :: i, j, k
    
        A(1:n,1:n) = Ain(1:n,1:n)
        L=0.d0
        U=0.d0
        b=0.d0
        
        do k=1, n-1
            do i=k+1,n
                coeff=a(i,k)/a(k,k)
                L(i,k) = coeff
                A(i,k+1:n) = A(i,k+1:n)-coeff*A(k,k+1:n)
            end do
        end do
    
        forall (i=1:n)  L(i,i) = 1.d0
        forall (j=1:n) U(1:j,j) = A(1:j,j)
        do k=1,n
            b(k)=1.d0
            d(1) = b(1)
            do i=2,n
                d(i)=b(i)
                d(i) = d(i) - dot_product(L(i,1:i-1),d(1:i-1))
            end do
            x(n)=d(n)/U(n,n)
            do i = n-1,1,-1
                x(i) = d(i)
                x(i)=x(i)-dot_product(U(i,i+1:n),x(i+1:n))
                x(i) = x(i)/U(i,i)
            end do
            A_inverse(1:n,k) = x(1:n)
            b(k)=0.d0
        end do
    
      end subroutine abq_inverse_LU
    