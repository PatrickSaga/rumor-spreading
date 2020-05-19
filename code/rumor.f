	program rumor_spreading

	
	! Rumor spreading model (SIHR) 


	! The model consists on a set of particles in one of the 4 possible states : 
	! I - ignorants
	! S - spreaders
	! H - hiberantors
	! R - stiflers

	! Particles evolve according to a master equation obtained under mean-field approximation.
	! The program runs a MC method using a Gillespie algorithm

	
	implicit none

	integer N								! total number of particles
	double precision ign, spr, hib, sti 				! fraction of particles in each state
	double precision t, tmax, t_trans, t_prev, dt 			! time variables and time-step between measures
	double precision k							! average degree of the network
	double precision u, dran_u 						! random number generator
	double precision lambda, beta, delta, xi, eta, alpha		! transition probabilities
	double precision rate1,rate2,rate3,rate4				! transition rates
	double precision rate5,rate6,rate7,rate8,rate			! and more transition rates
	

	! Parameters
	!lambda = 0.8d0
	!beta = 0.2d0
	!delta = 0.5d0
	!xi = 0.0d0
	!eta = 0.5d0
	!alpha = 0.5d0

	! Nekovee
	lambda = 0.09d0
	beta = 0.0d0
	delta = 0.9d0
	xi = 0.5d0
	eta = 0.0d0
	alpha = 0.5d0

	! SIR
	!lambda = 0.6d0
	!beta = 0.0d0
	!delta = 0.0d0
	!xi = 0.0d0
	!eta = 0.0d0
	!alpha = 0.5d0

	! Homeopathy
	!lambda = 5.5d0
	!beta = 1d0
	!delta = 0.5d0
	!xi = 7d0
	!eta = 9d0
	!alpha = 0.5d0

	t = 0.0d0
	tmax = 10000000.0d0
	dt = 100.0d0
	k = 100d0
	N = 1000000

	! Initial conditions
	ign = dble(N-1)/N
	spr = 1.0d0/N
	hib = 0.0d0
	sti = 0.0d0
	call dran_ini(481992)

	! open the output files
	open(unit=21,file='crit_thresh',status='replace',action='write')
	! write(21,*) "time      ignorants      spreaders      hibernators	stiflers"
	write(21,*) t, ",", ign,",", spr,",", hib,",", sti

	! Gillespie algorithm
	t_prev = 0.d0
	do while (t .lt. tmax)
	  rate1 = lambda*k*ign*spr
	  rate2 = beta*k*ign*spr
	  rate3 = delta*spr
	  rate4 = xi*hib
	  rate5 = eta*k*spr*hib
	  rate6 = alpha*k*spr*spr
	  rate7 = alpha*k*spr*hib
	  rate8 = alpha*k*spr*sti
	  rate = rate1+rate2+rate3+rate4+rate5+rate6+rate7+rate8

	  t_trans = -dlog(dran_u())/rate
	  t = t+t_trans
	  u = dran_u()*rate

	  if (u < rate1) then
	    ign = ign-1.d0/N
	    spr = spr+1.d0/N
	  elseif (u < rate1+rate2) then
	    ign = ign-1.d0/N
	    sti = sti+1.d0/N
	  elseif (u < rate1+rate2+rate3) then
	    spr = spr-1.d0/N
	    hib = hib+1.d0/N
	  elseif (u < rate1+rate2+rate3+rate4+rate5) then
	    hib = hib-1.d0/N
	    spr = spr+1.d0/N
	  else
	    spr = spr-1.d0/N
	    sti = sti+1.d0/N
	  endif
	  if (t-t_prev > dt) then
	    write(21,*) t, ",", ign,",", spr,",", hib,",", sti
	    t_prev = t
	  endif

	enddo
	close(21)
	end program rumor_spreading
