implicit none

integer,parameter:: step_number=100000, grid=50, trajectory=2000
real:: xmin,xmax,pmin,pmax,m,a,b,c,d,dp,p,vel,pstart,x,xt,velt,v11,v12,v22
real:: dt,t,ttot,qcl,pcl,k_mn,k_mx
real:: vnew,vcor,d_ij,ener
real:: dve,nonad_coup,random,hopping,transm,reflect
real:: p0(grid+1),ti(step_number),adiab_ener1,adiab_ener2
real:: a_11,a_12,a_21,a_22,b_12,b_21,p_12,p_21
integer:: adiab_surface,i,j,k,l

complex:: ci_0(2),ci_t(2),ct

open(3000,file='data.out')



xmin = -10
xmax = 10
pmin = 0
pmax = 35

m = 2000
a = 0.01
b = 1.6
c = 0.005
d = 1

dp = (pmax-pmin)/grid
do i = 1, grid+1
p0(i) = (i-1)*dp
end do

dt = 1.0      !timestep
do j = 1, step_number
ti(j) = j*dt
end do

!----------------------------------------

do k = 1, grid

pstart = p0(k)
hopping=0
transm = 0
reflect = 0

do j = 1, trajectory
!write(*,*) j  
adiab_surface=1
x = -10.d0
p = pstart
vel = p/m
ci_t(1) = (1.d0,0.d0)
ci_t(2) = (0.d0,0.d0)

!-------------------------------------

do i = 1, step_number

t = ti(i)

call nonadiabatic_coupling(x,p,nonad_coup)

d_ij = nonad_coup
!write(101,*) t, d_ij
do l = 1, 2
call rk4(ci_t(1),ci_t(2),x,vel,dt,ct,l,d_ij)
ci_t(l) = ct
end do
!write(102,*) t, ci_t(1),ci_t(2)

call vel_ver(adiab_surface,t,dt,x,p,xt,velt)

x = xt
p = m*velt

a_11= cabs(ci_t(1))
a_22= cabs(ci_t(2))
a_12= ci_t(1)*conjg(ci_t(2))
a_21= ci_t(2)*conjg(ci_t(1))
b_21= 2*real(d_ij*a_12*velt)
b_12= -2*real(d_ij*a_21*velt)
p_21= (b_21/a_11)*dt
p_12= (b_12/a_22)*dt

!write(100,*) t, b_21,b_12,p_21,p_12

if(p_21<0)then
p_21 = 0
end if

if(p_12<0)then
p_12 = 0
end if

call random_number(random)


        if (adiab_surface==1) then

         if(p_21>random)then
           dve = adiab_ener2(x) - adiab_ener1(x)

              if (dve>(p**2/(2*m))) then
                adiab_surface = 1
              else
                adiab_surface = 2
                vcor = sqrt((2/m)*((p**2/(2*m))+(-dve)))

                  if (p<0) then
                   vnew = -vcor
                  else
                   vnew = vcor
                  end if

                p = m*vnew
              end if

         else
          adiab_surface = 1
         end if

        else if (adiab_surface==2) then

         if (p_12>random) then

           adiab_surface = 1
           dve = adiab_ener2(x) - adiab_ener1(x)
           vcor = sqrt((2/m)*(dve + (p**2/(2*m))))

              if (p<0) then
                vnew = -vcor
              else
                vnew = vcor
              end if

             p = m*vnew

         else
          adiab_surface = 2
         end if

        end if

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (x<xmin) then
 if (adiab_surface==1) then
  reflect = reflect + 1
 else
  hopping = hopping + 1
 end if
 exit
end if
if (x>xmax) then
 if (adiab_surface==1) then
  transm = transm + 1
 else
  hopping = hopping + 1
 end if
 exit
end if

!-----------------------------------------------------------------------

end do

end do

write(3000,*) pstart, transm, reflect, hopping

end do

end

!-----------------------------------------------------------------------

FUNCTION V11(x)
implicit none
real:: m,a,b,c,d,v11,x
m = 2000
a = 0.01
b = 1.6
c = 0.005
d = 1

if (x<0) then
v11 = -a*(1-(exp(b*x)))
else
v11 = a*(1-(exp(-b*x)))
end if

END FUNCTION v11

!---------------------------------------

FUNCTION V22(x)
implicit none
real:: m,a,b,c,d,x,v22,v11
m = 2000
a = 0.01
b = 1.6
c = 0.005
d = 1
v22 = -v11(x)

END FUNCTION v22

!----------------------------------------------------------

FUNCTION V12(x)
implicit none
real:: m,a,b,c,d,x,v12
m = 2000
a = 0.01
b = 1.6
c = 0.005
d = 1
v12 = c*(exp(-d*(x**2)))

END FUNCTION v12

!-------------------------------------------------------------

FUNCTION adiab_ener1(x)
implicit none
real:: adiab_ener1,x,trace,det,v11,v22,v12

trace = 0
det = (v11(x)*v22(x)) - (v12(x)*v12(x))
adiab_ener1 = -sqrt(-(4*det))/2

END FUNCTION adiab_ener1

!------------------------------------------------------------------

FUNCTION adiab_ener2(x)
implicit none
real:: adiab_ener2,x,trace,det,v11,v22,v12

trace = 0
det = (v11(x)*v22(x)) - (v12(x)*v12(x))
adiab_ener2 = sqrt(-(4*det))/2

END FUNCTION adiab_ener2


!               SUBROUTINE NONADIABATIC COUPLING

subroutine nonadiabatic_coupling(x1,p,d12)
implicit none
real:: a,b,c,d,m,v11,v22,v12,dv11,dv22,dv12,x1,p,d12
real:: c11,c22,c12,c21,dx1f11,dx1f22,dx1g,d12_x1
integer,parameter:: dimn = 2, LWORK=(3*dimn + 1), LDA=dimn

real(8):: ham(dimn,dimn),w(dimn),work(LWORK),b1(dimn,dimn),s(dimn,dimn),c1(dimn)
integer:: info

m = 2000
a = 0.01
b = 1.6
c = 0.005
d = 1
if (x1.lt.0) then
v11 = -a*(1-(exp(b*x1)))
v22 = -v11
dv11 = a*b*exp(b*x1)
dv22 = -dv11
else
v11 = a*(1-(exp(-b*x1)))
v22 = -v11
dv11 = a*b*exp(-b*x1)
dv22 = -dv11
end if

v12 = c*exp(-d*(x1**2))
dv12 = -2*d*x1*v12

ham(1,1) = v11
ham(2,2) = v22
ham(1,2) = v12
ham(2,1) = v12

dx1f11 = dv11
dx1f22 = dv22
dx1g = dv12

 call DSYEV('V','L',dimn,ham,LDA,W,WORK,LWORK,info)

c11 = ham(1,1)
c21 = ham(2,1)
c12 = ham(1,2)
c22 = ham(2,2)

d12_x1 = c11*c12*dx1f11 + c11*c22*dx1g + c21*c12*dx1g + c21*c22*dx1f22

d12 = ((d12_x1))/(w(2)-w(1))

end subroutine nonadiabatic_coupling

!velocity verlet

subroutine vel_ver(surface,t,dt,x,p,xt,vt)
implicit none
real:: dt,t,x,p,xt,vt,acc,acc1,m,a,b,c,d,f_p
integer:: surface,i,j

m = 2000
a = 0.01
b = 1.6
c = 0.005
d = 1
acc = f_p(t,x,p,surface)/m

xt = x + (p/m)*dt + 0.5*acc*(dt**2)

acc1 = f_p(t+dt,xt,p,surface)/m

vt = (p/m) + 0.5*(acc+acc1)*dt

end subroutine vel_ver

!--------------------------------------------

function f_p(t,r,p,surf)
implicit none
real:: r, p, f_p,t,a,b,c,d,m,mag_f,v11,v22,v12,v11p,v22p,v12p
integer:: surf


m = 2000
a = 0.01
b = 1.6
c = 0.005
d = 1

if (r.lt.0) then
 v11 = -a*(1-exp(b*r))
 v11p = a*b*exp(b*r)
 v12 = c*exp(-d*(r**2))
 v12p = -2*d*r*v12
else
 v11 = a*(1-exp(-b*r))
 v11p = a*b*exp(-b*r)
 v12 = c*exp(-d*(r**2))
 v12p = -2*d*r*v12
end if

mag_f = ((v11*v11p)+(v12*v12p))/sqrt((v11**2)+(v12**2))

if (surf.eq.1) then
f_p = mag_f
else if (surf.eq.2) then
f_p = -mag_f
end if

end function f_p


!runge-kutta

subroutine rk4(c1,c2,x,v,dt,ct,neq,d12)
implicit none

complex,parameter:: iota=(0,1)
complex:: c1,c2,ct,f_c1,f_c2,k1,k2,k3,k4
integer:: neq
real:: x,v,p,a,b,c,d,m,d12,dt,t

m = 2000
a = 0.01
b = 1.6
c = 0.005
d = 1
p = v*m

if (neq==1) then
 k1 = f_c1(t,c1,c2,x,p,d12)
 k2 = f_c1((t+dt/2),(c1+(dt*k1/2)),(c2+(dt*k1/2)),x,p,d12)
 k3 = f_c1((t+dt/2),(c1+(dt*k2/2)),(c2+(dt*k2/2)),x,p,d12)
 k4 = f_c1((t+dt),(c1+(dt*k3)),(c2+(dt*k3)),x,p,d12)

 ct = c1 + (dt/6)*(k1 + (2*k2) + (2*k3) + k4)
else
 k1 = f_c2(t,c1,c2,x,p,d12)
 k2 = f_c2((t+dt/2),(c1+(dt*k1/2)),(c2+(dt*k1/2)),x,p,d12)
 k3 = f_c2((t+dt/2),(c1+(dt*k2/2)),(c2+(dt*k2/2)),x,p,d12)
 k4 = f_c2((t+dt),(c1+(dt*k3)),(c2+(dt*k3)),x,p,d12)

 ct = c2 + (dt/6)*(k1 + (2*k2) + (2*k3) + k4)
end if

end subroutine rk4
!                         FUNCTION f_c1
function f_c1(t,c1,c2,r,p,d12)
implicit none

complex,parameter:: iota=(0,1)
complex:: c1,c2,f_c1
real::a,b,c,d,m,r,p,v,d12,t,adiab_ener1
m = 2000
a = 0.01
b = 1.6
c = 0.005
d = 1
v = p/m

f_c1 = -iota*c1*adiab_ener1(r) - c2*v*d12

end function f_c1
!                         FUNCTION f_c2
function f_c2(t,c1,c2,r,p,d12)
implicit none

complex,parameter:: iota=(0,1)
complex:: c1,c2,f_c2
real::a,b,c,d,m,r,p,v,d12,t,adiab_ener2
m = 2000
a = 0.01
b = 1.6
c = 0.005
d = 1
v = p/m

f_c2 = -iota*c2*adiab_ener2(r) + c1*v*d12

end function f_c2



