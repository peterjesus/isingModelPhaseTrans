*   CONDENSED MATTER THEORY
*   ASSIGNMENT 3. EXERCISE 2
*   
*   PHASE TRANSITION IN THE 2D-ISING MODEL
*   2Dising.f
************************************************************************
*
*   Script for calculate several quantities related with the
*   2-dimensional Ising model
*
*     a) Calculate E/J and M for each spin configuration
*     b) For each M, calculate Z and F(M,T)/J 
*     c) Plot F(M) for each value of T
*
************************************************************************
      
      
      parameter(n=4)              !Size of our lattice
      parameter(nsite=n*n)        !Number of sites in our lattice
      parameter(nconf=2**nsite)   !Number of possible configurations
      real*8 o,en,mag,v,e,m
      real*8 T,h,w,z,f
      dimension o(0:nsite-1),v(0:nsite-1)   !0 is the upper-left site and nsite-1 is the lower-right side
      dimension e(0:nconf-1),m(0:nconf-1)

!a) Calculate E/J and M for each spin configuration
      open(10,file='2Dising-A.dat')
      write(10,*) 'Conf(Decimal)       ','Energy/J           ','Magnetiz
     &ation'
     
*-initialization of vectors that will contains the values of energy/J and magnetization
      do l=0,nconf-1
         e(l)=0.0d0
         m(l)=0.0d0
      end do

*-main loop
      do l=0,nconf-1

*---initialization of vectors v and o for our spin-conf
        do i=0,nsite-1
           o(i)=0.0d0
           v(i)=0.0d0
        end do
      
*---going from decimal to binary
        ll=l
        do k=0,nsite-1              !looking for all the values smaller than the number of sites
           kk=nsite-1-k
           if((ll-2**kk).ge.0) then
             ll=ll-2**kk               !we rest the potence of 2 considered
             v(abs(nsite-1-kk))=1.0d0  !we memorize in our vector this point
             if(ll.eq.0) then
               go to 100               !if our final number is zero: work done
             end if
           else
             continue
           end if
        end do
        
*---getting the vector of spins for our configuration (all +1 or -1)
100     do i=0,nsite-1           !lattice
            o(i)=2.0d0*v(i)-1.0d0
        end do

*---initialization of energy and magnetization for each configuration
        en=0.0d0
        mag=0.0d0

*---calculate of the energy and the magnetization for our configuration
        do i=0,n-1
          do j=0,n-1
            if(j.eq.n-1 .and. i.eq.n-1) then
               en=en-o(i+n*j)*o(i)-o(i+n*j)*o(n*j)
            else
              if(j.eq.n-1) then
               en=en-o(i+n*j)*o(i+n*j+1)-o(i+n*j)*o(i)
              elseif(i.eq.n-1) then
               en=en-o(i+n*j)*o(n*j)-o(i+n*j)*o(i+n*(j+1))
              else
               en=en-o(i+n*j)*o(i+n*j+1)-o(i+n*j)*o(i+n*(j+1))
              endif
            endif
            mag=mag+o(i+n*j)
          end do
        end do
        e(l)=en
        m(l)=mag
        
      end do
      
      do l=0,nconf-1
         write(10,*) l,e(l),m(l)
      end do
      
      close(10)
      
      
!b)For each M, calculate Z and F(M,T)/J
      open(20,file='2Dising-B-3p5.dat')
      open(30,file='2Dising-B-3p6.dat')
      open(40,file='2Dising-B-3p7.dat')
      open(50,file='2Dising-B-3p8.dat')
      
*-definition of the possible values of kT/J
      ntemp=3
      T=3.5d0
      h=0.1d0
      
*-main loop
      do nn=0,ntemp
      
        T=T+h*nn 
        mw=10*(2+nn)
        
        write(mw,*) 'Magnetization  ','Partition Function          ','Fr
     &ee Energy'
              
        do nmag=-nsite,nsite,2      !Possible values of Magnetiz: min,max val; even val (in order to avoid non-possible cases)
        
*---initialization of partition function and the free energy
          z=0.0d0
          f=0.0d0
        
*---now we want to look for a configuration has the same magnetiz-value that the value
          do l=0,nconf-1
          
            if(m(l).eq.nmag) then     !if the condition is satisfied, we take the number l associated to the spin conf
            
*-------initialization of vectors v and o for our spin-conf
              do i=0,nsite-1
                o(i)=0.0d0
                v(i)=0.0d0
              end do
            
*-------going from decimal to binary
              ll=l
              do k=0,nsite-1              !looking for all the values smaller than the number of sites
                  kk=nsite-1-k
                if((ll-2**kk).ge.0) then
                  ll=ll-2**kk               !we rest the potence of 2 considered
                  v(abs(nsite-1-kk))=1.0d0  !we memorize in our vector this point
                  if(ll.eq.0) then
                    go to 200               !if our final number is zero: work done
                  end if
                else
                  continue
                end if
              end do

*-------getting the vector of spins for our configuration (all +1 or -1)
200           do i=0,nsite-1           !lattice
                o(i)=2.0d0*v(i)-1.0d0
              end do
              
*-------we know that the partition function is the exponential of -beta*hamiltonian evaluated for our configuration
*-------we know that our partition function is the sume (for all the states with an magnetization M)
*-------of the exponential of -beta*hamiltonanian (where the hamiltonian -or energy- is evaluated for this states)
              do i=1,n
                do j=0,n-1
                  if(j.eq.n-1 .and. i.eq.n-1) then
                        w=o(i+n*j)*o(i)-o(i+n*j)*o(n*j)
                  else
                    if(j.eq.n-1) then
                        w=o(i+n*j)*o(i+n*j+1)-o(i+n*j)*o(i)
                    elseif(i.eq.n-1) then
                        w=o(i+n*j)*o(n*j)-o(i+n*j)*o(i+n*(j+1))
                    else
                        w=o(i+n*j)*o(i+n*j+1)-o(i+n*j)*o(i+n*(j+1))
                    endif
                  endif
                  z=z+exp(-w/T)     !after get the energy for the state, we sum the exponential
                end do
              end do
            
*-------the free energy is -beta^(-1)*ln(partition function)
              f=-T*log(z)

            else
              continue
            endif
            
          end do
        
*---now we write the partition function and the free energy obtained for a certain magnetization
          write(mw,*) nmag,z,f
          
          
          
        end do
      end do
      
      close(20)
      close(30)
      close(40)
      close(50)
      
      stop
      end
