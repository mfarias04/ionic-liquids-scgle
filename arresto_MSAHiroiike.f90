!dinamica de la SCGLE para calcular el desplazamiento cuadratico medio y las D0 de dos especies ionicas
!
!
program kasuo_h
implicit none
Integer, Parameter      :: nesp= 2, nkmax = 2**16
Real(8), Parameter      :: pi= dacos(-1.d0), k0 = .0001d0, dk = 0.001d0
Real(8), Parameter	:: fact =1.d0/(6.d0*(pi**2))
integer			:: nk, l, j , it, i, m
Real(8), Dimension(nesp,nesp,nkmax)	:: lamdaI, sqrtnh, csqrtn, Q, QI, lambdaI
Real(8)			:: k, suma1, suma2, suma3, tol, error, w, D, gama, suma4, kvw
Real(8)			:: a(nesp), b(nesp), psi(3), Ts, tolH, large, Tmax, Tmin, phistar, phivw
Real(8), Dimension(nesp)   :: X, rho, z, sigma, N, kc, errorH, suma, errorFH, rhovw
Real(8), Dimension(nesp,nesp) :: C2, C1, C, F, F1, F2, F3, F4, IM, Rhoa, RhoaI, Cbaxter, cbax, gm, Rhoavw, RhoaIvw
Real(8), Dimension(nesp,nesp) :: chiroike, ctotal, S, SInv, Sh, ShI, one, oneI, two, twoI, three, Id, gamma
Real(8), Dimension(nesp, 10000) :: gmlist
integer, dimension(nesp)	::flag
double precision, Dimension(nesp) :: flag1

open(unit = 1,file = 'sk_hiro.dat', status = 'unknown')

z=(/1.d0, -1.d0/)			!carga de las partículas
sigma=(/1.d0,0.1d0/)			!diametro de las partículas


IM= 0.0d0
do l=1, nesp

   IM(l,l)=1.0d0

end do

!Ts=1.d4
write(*,*) "Fraccion de volumen = "
read(*,*) phistar


write(*,*) "Temperatura reducida = "
read(*,*) Ts

!phistar=0.2


rho(1)= -6.d0*phistar*z(2)/(pi*(z(1)*(sigma(2))**3 - z(2)*(sigma(1))**3))

rho(2)= 6.d0*phistar*z(1)/(pi*(z(1)*(sigma(2))**3 - z(2)*(sigma(1))**3))


Rhoa=0.d0
RhoaI=0.d0

do l=1,nesp
do j=1,nesp

Rhoa(l,l)=sqrt(rho(l))
RhoaI(l,l)=1.d0/sqrt(rho(l))

enddo
enddo


call sdkhiroike(Q,QI, lambdaI)


do i = 1, nkmax

sqrtnh(:,:,i) = matmul(Q(:,:,i) - IM, RhoaI)
csqrtn(:,:,i) = matmul(RhoaI, IM - QI(:,:,i))

!	write(*,*) SkI(1,1,1)

end do


call arresto(2,nkmax,rho,lambdaI,Q,QI,sqrtnh,csqrtn,gamma,flag1)


close(unit = 1)


contains


subroutine sdkhiroike(Q,QI,lambdaI)	!!!!calculando la matriz de factores de estructura de hiroike + baxter
  double precision :: Q,QI
  double precision :: lambdaI
integer	:: nk, it, m
real(8)	:: k, kc

dimension Q(2,2,nkmax),QI(2,2,nkmax), lambdaI(2,2,nkmax)


lamdaI=0.0d0
lambdaI=0.d0
kc = 2*pi*1.305d0!3.2d0!2*pi*1.305d0!0.939d0
do nk = 1, nkmax

	k = k0 + (nk-1)*dk

call baxter(k, S)


call inversion(S, SInv)

	Cbaxter= IM - SInv

call Correlacion(k, C)

	cbax=matmul(matmul(RhoaI,Cbaxter),RhoaI)
	chiroike= C + cbax
	ctotal=matmul(matmul(Rhoa,chiroike),Rhoa)

	ShI(:,:)= IM - ctotal

call inversion(ShI(:,:), Sh(:,:))

!!!!
! Sh=S
!!!!
 Q(:,:,nk)= Sh


call inversion(Q(:,:,nk), QI(:,:,nk))

do m = 1, 2
   lambdaI(m,m,nk) = (1.d0 + (k/kc)**2)
end do

	write(1,*) sngl(k), sngl(Q(1,1,nk)),sngl(Q(2,2,nk))

end do




end subroutine sdkhiroike

subroutine gama_hiroike			!aqui se encuentra la gama de hiroike por iteración


suma1=0.0d0
do j=1, nesp
	suma1=suma1 + rho(j)*(sigma(j))**3
end do

w= (pi/2.d0)*(1.d0/(1.d0 - (pi/6.d0)*suma1))			!aqui se obtiene la "c"


gama=1.d-8
tol=1.d-6
error=1.d0
it=0


do while ( error > tol)

	it= it + 1

	suma2=0.0d0
	suma3=0.0d0
	do j=1, nesp

		suma2=suma2 + rho(j)*sigma(j)*z(j)*(1.d0/(1.d0 + gama*sigma(j)))
		suma3=suma3 + rho(j)*((sigma(j))**3)*(1.d0/(1.d0 + gama*sigma(j)))
	end do


do l=1,nesp

X(l)= z(l)/(1.d0 + gama*sigma(l)) - (w*(sigma(l))**2)*suma2/((1.d0 + w*suma3)*(1.d0 + gama*sigma(l)))

end do



	do l=1 , nesp

		N(l) = (X(l) - z(l))/sigma(l)
	end do



	suma4=0.0d0
	do l=1, nesp
		suma4= suma4 + rho(l)*(X(l))**2
	end do

D= suma4


	error= abs(((sqrt(pi/Ts))*(sqrt(D)) - gama)/gama)

	gama=(sqrt(pi/Ts))*(sqrt(D))

!	print*, "it, gama=", it, sngl(gama)
end do


end subroutine gama_hiroike



subroutine Correlacion(k, C)		!calculando la funcion de correlacion de hiroike

Real(8), Dimension(nesp,nesp) :: C, F, F1, F2, F3, F4
integer		:: l,j
Real(8)		:: k



call gama_hiroike


do l=1, nesp

	do j=1, nesp
F(l,j)= -z(l)*N(j) + X(l)*(N(l) + gama*X(l)) - (sigma(l)/3.d0)*(N(l) + gama*X(l))**2

F1(l,j)=(sigma(l) - sigma(j))*(((X(l) + X(j))/4.d0)*(N(l) + gama*X(l) - N(j) - gama*X(j)) -&
((sigma(l) -sigma(j))/16.d0)*((N(l) + gama*X(l) + N(j) + gama*X(j))**2 - 4.d0*N(l)*N(j)))

F2(l,j)=(X(l) - X(j))*(N(l) - N(j)) + ((X(l))**2 + (X(j))**2)*gama + (sigma(l) + sigma(j))*&
N(l)*N(j) - (1.d0/3.d0)*(sigma(l)*(N(l) + gama*X(l))**2 + sigma(j)*(N(j) + gama*X(j))**2)

F3(l,j)= (X(l)/sigma(l))*(N(l) + gama*X(l)) + (X(j)/sigma(j))*(N(j) + gama*X(j)) +&
 N(l)*N(j) - (0.5d0)*((N(l) + gama*X(l))**2 + (N(j) + gama*X(j))**2)

F4(l,j)=((N(l) + gama*X(l))**2)/(6.d0*(sigma(l)**2)) + ((N(j) + gama*X(j))**2)/(6.d0*&
(sigma(j)**2))



	end do

end do


do l=1, nesp
	do j=1, nesp

if(l==j) then
	C1(l,j)=0.0d0
Else If (sigma(l) < sigma(j)) Then
	C1(l,j)=(4.d0*F(l,j)/((k**2)*Ts*sqrt(2.d0*pi)))*((abs((sigma(l) - sigma(j))/2.d0))*&
cos(k*abs((sigma(l) - sigma(j))/2.d0)) - (1.d0/k)*sin(k*abs((sigma(l) - sigma(j))/2.d0)))
end if


 C2(l,j)=(2.d0*F1(l,j)/((k**2)*Ts*sqrt(2.d0*pi)))*(cos(k*abs((sigma(l) - sigma(j))/2.d0)) -&
cos(k*(sigma(l) + sigma(j))/2.d0)) - (2.d0*F2(l,j)/(k*Ts*sqrt(2.d0*pi)))*&
((1.d0/k)*(abs((sigma(l)- sigma(j))/2)*cos(k*abs((sigma(l) - sigma(j))/2)) - ((sigma(l) +&
sigma(j))/2.d0)*cos(k*(sigma(l) + sigma(j))/2.d0)) + (1.d0/k**2)*(sin(k*((sigma(l) +&
sigma(j))/2.d0)) - sin(k*abs((sigma(l) - sigma(j))/2)))) + (2.d0*F3(l,j)/(Ts*k*sqrt(2.d0*pi)))*&
((1.d0/k)*(((abs((sigma(l) - sigma(j))/2.d0))**2)*cos(k*abs((sigma(l) - sigma(j))/2.d0)) -&
(((sigma(l) + sigma(j))/2.d0)**2)*cos(k*(sigma(l) + sigma(j))/2.d0)) + (2.d0/k**2)*(((sigma(l) +&
 sigma(j))/2.d0)*sin(k*(sigma(l) + sigma(j))/2.d0) - (abs((sigma(l) - sigma(j))/2.d0))*sin(k*&
abs((sigma(l) - sigma(j))/2.d0))) + (2.d0/k**3)*(cos(k*(sigma(l) + sigma(j))/2.d0) -&
 cos(k*abs((sigma(l) - sigma(j))/2.d0)))) + (2.d0*F4(l,j)/(k*Ts*sqrt(2.d0*pi)))*((1.d0/k)*&
(((abs((sigma(l) - sigma(j))/2.d0))**4)*cos(k*abs((sigma(l) - sigma(j))/2.d0)) - (((sigma(l) +&
sigma(j))/2.d0)**4)*cos(k*(sigma(l) + sigma(j))/2.d0)) + (4.d0/k**2)*((((sigma(l) +&
 sigma(j))/2.d0)**3)*sin(k*(sigma(l) + sigma(j))/2.d0) - ((abs((sigma(l) - sigma(j))/2.d0))**3)*&
sin(k*abs((sigma(l) - sigma(j))/2.d0))) + (12.d0/k**3)*((((sigma(l) + sigma(j))/2.d0)**2)*&
cos(k*(sigma(l) + sigma(j))/2.d0) - ((abs((sigma(l) - sigma(j))/2.d0))**2)*cos(k*abs((sigma(l) -&
 sigma(j))/2.d0))) + (24.d0/k**4)*((abs((sigma(l) - sigma(j))/2.d0))*sin(k*&
abs((sigma(l) - sigma(j))/2.d0)) - ((sigma(l) + sigma(j))/2.d0)*sin(k*(sigma(l) +&
 sigma(j))/2.d0)) + (24.d0/k**5)*(cos(k*abs((sigma(l) - sigma(j))/2.d0)) - cos(k*(sigma(l) +&
 sigma(j))/2.d0))) - (2.d0*z(l)*z(j)/((k**2)*Ts*sqrt(2.d0*pi)))*cos(k*(sigma(l) +&
 sigma(j))/2.d0)


	end do
end do

Do l = 1,nesp
	Do j = 1,nesp

		If (sigma(l) > sigma(j)) Then
			C1(l,j) = C1(j,l)
		End If
	End Do
End Do

	C = ((sqrt(2.d0*pi))**3)*(C1 + C2)


end subroutine Correlacion


subroutine inversion(a,b)
	Real(8)			:: det
	Real(8), dimension(nesp,nesp), intent(in)	:: a
	Real(8), dimension(nesp,nesp), intent(out)	:: b

	det = a(1,1)*a(2,2)-a(1,2)*a(2,1)

	b(1,1) = a(2,2)/det
	b(1,2) = -a(1,2)/det
	b(2,1) = -a(2,1)/det
	b(2,2) = a(1,1)/det
return

end subroutine


subroutine baxter(k, S)		!calculando el factor de estructura para mezclas de baxter
complex(8), dimension(nesp,nesp) :: QM, QMT
Real(8), Dimension(nesp, nesp)	:: S, SInv
Complex(8), Parameter	:: I = Cmplx(0.d0,1.d0)
integer 		:: l,j
Real(8)			:: k




do j=1, 3

	psi(j)= (pi/6.d0)*dot_product(rho,sigma**j)

end do

do  l=1, nesp

	a(l)= (1.d0 - psi(3) + 3.d0*sigma(l)*psi(2))/(1.d0 - psi(3))**2

	b(l)=-((3.d0/2.d0)*((sigma(l))**2)*psi(2))/(1.d0 - psi(3))**2

end do




do l=1, nesp
	do j=1,nesp

	QM(l,j)= IM(l,j) - (pi*(sqrt(rho(l)*rho(j)))/k**3)*(2.d0*exp(I*k*((sigma(l) +&
sigma(j))/2.d0))*(k*(b(l) + a(l)*((sigma(l) + sigma(j))/2.d0)) + a(l)*I) + exp(I*k*((sigma(l)-&
 sigma(j))/2.d0))*(a(l)*(k**2)*I*(((sigma(l)-sigma(j))/2.d0)**2 - ((sigma(j) +&
 sigma(l))/2.d0)**2)- 2.d0*I*a(l) - 2.d0*b(l)*I*(k**2)*(sigma(j)) -2.d0*k*(a(l)*&
((sigma(l)-sigma(j))/2.d0) + b(l))))

	QMT(l,j)= IM(l,j) + (pi*(sqrt(rho(l)*rho(j)))/k**3)*(2.d0*exp(I*(-k)*((sigma(j) +&
sigma(l))/2.d0))*((-k)*(b(j) + a(j)*((sigma(j) + sigma(l))/2.d0)) + a(j)*I) + exp(I*(-k)*&
((sigma(j)- sigma(l))/2.d0))*(a(j)*(k**2)*I*(((sigma(j)-sigma(l))/2.d0)**2 - ((sigma(l)+&
 sigma(j))/2.d0)**2) - 2.d0*I*a(j) - 2.d0*b(j)*I*(k**2)*(sigma(l)) - 2.d0*(-k)*(a(j)*&
((sigma(j)-sigma(l))/2.d0) + b(j))))
	end do
end do


SInv = matmul(QMT,QM)


call inversion(SInv, S)

return

end subroutine baxter



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Subrutina de arresto de Ernie

subroutine arresto(nm,nk,rho,lambdaI,Sk,SkI,sqrtnh,csqrtn,gamma,flag)
integer	:: nk, it, i, nm
double precision :: k,suma,func,func1,func2,func3,fact
double precision :: error,tol,max,SkI,gammai,gamma,Id,Sk
double precision :: lambda,lambdaI,rho,rhoa,rhoaI,flag
double precision :: func1I,func2I,func3I,gammalist,sqrtnh,csqrtn


dimension gamma(nm,nm),error(nm),Id(nm,nm),func(nm)
dimension rho(nm),rhoa(nm,nm),rhoaI(nm,nm),suma(nm)
dimension func1(nm,nm),func1I(nm,nm),func2(nm,nm),func2I(nm,nm)
dimension func3(nm,nm),func3I(nm,nm),gammalist(nm,10000)
dimension lambdaI(nm,nm,nk),SkI(nm,nm,nk),flag(nm)
dimension sqrtnh(2,2,nkmax), csqrtn(2,2,nkmax), Sk(nm,nm,nk)

! write(*,*) Sk(1,1,1), lambdaI(1,1,1), sqrtnh(1,1,1), csqrtn(1,1,1), dk

func1 = 0.d0
func2 = 0.d0
func3 = 0.d0
func1I = 0.d0
func2I = 0.d0
func3I = 0.d0
flag=0.d0
gammalist=0.d0

rhoa = 0.d0
rhoaI = 0.d0

     do i = 1, nm
   	rhoa(i,i) = sqrt(rho(i))
   	rhoaI(i,i) = 1.d0/rhoa(i,i)
     end do
!write(*,*) rhoa

!definicion de la gama
gamma = 0.d0
     do l = 1, nm
  	gamma(l,l) = 1.d-6
     end do

Id = 0.d0
     do l = 1, nm
   	Id(l,l) = 1.d0
     end do

fact = 1.d0/(6.d0*(pi**2))
tol = 1.d-6
error = 1.d0
it = 0
max = 1.d+33


     do while(error(1)>tol .and. error(2)>tol .and. gamma(1,1)<max .and. gamma(2,2)<max)
	 it = it+1
	 suma = 0.d0
	 do i = 1, nk
		k = k0 + (i-1)*dk !dk*i
		func1 = Id + k**2*matmul(gamma,lambdaI(:,:,i))
		call INVERSION1(nm,func1,func1I)
		func2= Id + k**2*matmul(matmul(gamma,lambdaI(:,:,i)),SkI(:,:,i))
		call INVERSION1(nm,func2,func2I)
		func3=matmul(matmul(csqrtn(:,:,i),func2I),sqrtnh(:,:,i))
!                write(*,*) sqrtnh(:,:,i)
!                write(*,*) csqrtn(:,:,i)
		do l = 1, nm
			suma(l)= suma(l) + (k**4)*func1I(l,l)*func3(l,l)
		end do
!                write(*,*) suma

	end do
        func = fact*suma*dk
	do l = 1, nm
		error(l) = abs((1.d0/func(l) - gamma(l,l))/gamma(l,l))
	end do

	do l = 1, nm
		gamma(l,l) = 1.d0/func(l)
	end do

	do l = 1, nm
		gammalist(l,it) = gamma(l,l)
	end do
        print*,"Iteration=",it ,"Gamma_i=",gamma(1,1), gamma(2,2)

     end do

	do l = 1, nm
		flag(l)= sign(gammalist(l,it) - 2.0d0*gammalist(l,it - 1) + gammalist(l,it - 2),1.d0)
!                flag(l)= gammalist(l,it) - 2.0d0*gammalist(l,it - 1) + gammalist(l,it - 2)
	end do

!write(*,*) "gamma(1,1) / gamma (2,2)"

! write(*,*) gamma(1,1), gamma(2,2)

end subroutine arresto




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subrutina de inversión de matrices

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! SE INCORPORA LA SUBROUTINE DE INVERSION QUE NOS PROPORCIONO ALEX
!---------------------------------------------------------------------------
!DI: NOS DETERMINA LA DIMENSION DE LA MATRIZ
!A:  ES LA MATRIZ ENTRANTE PARA QUE SE PUEDA INVERTIR
!B:  ES LA MATRIZ INVERSA
	SUBROUTINE INVERSION1(DI,A,B)
        implicit none
	INTEGER DI,INDX(DI), I, J, N
	REAL(8) A(DI,DI),B(DI,DI)
	REAL	AA(DI,DI), D, Y(DI, DI), BB(DI, DI)

        N=DI
	AA = A

	!       SET UP IDENTITY MATRIX.
	DO I=1,N
	DO J=1,N
	Y(I,J)=0.
	ENDDO
	Y(I,I)=1.
	ENDDO
	!                                        DECOMPOSE THE MATRIX JUST ONCE.
	CALL LUDCMP(AA,N,DI,INDX,D)
	!                                        FIND INVERSE BY COLUMNS.
	DO J=1,N
	CALL LUBKSB(AA,N,DI,INDX,Y(1,J))
	!        NOTE THAT FORTRAN STORES TWO-DIMENSIONAL MATRICES BY COLUMN, SO Y(1,J) IS THE
	!        ADDRESS OF THE JTH COLUMN OF Y.
	ENDDO
	B=Y
	RETURN
	END SUBROUTINE INVERSION1
	!

    !.................................................................

	!

	SUBROUTINE LUBKSB(A,N,DI,INDX,B)
        implicit none
	INTEGER N,DI,INDX(N)
	REAL A(DI,DI),B(N)
	!     SOLVES THE SET OF N LINEAR EQUATIONS A · X = B. HERE A IS INPUT, NOT AS THE MATRIX A BUT
	!     RATHER AS ITS LU DECOMPOSITION, DETERMINED BY THE ROUTINE LUDCMP. INDX IS INPUT AS THE
	!     PERMUTATION VECTOR RETURNED BY LUDCMP. B(1:N) IS INPUT AS THE RIGHT-HAND SIDE VECTOR B,
	!     AND RETURNS WITH THE SOLUTION VECTOR X. A, N, DI, AND INDX ARE NOT MODIϬED BY THIS ROUTINE
	!     AND CAN BE LEFT IN PLACE FOR SUCCESSIVE CALLS WITH DIϬERENT RIGHT-HAND SIDES B. THIS ROUTINE
	!     TAKES INTO ACCOUNT THE POSSIBILITY THAT B WILL BEGIN WITH MANY ZERO ELEMENTS, SO IT IS EϬCIENT
	!     FOR USE IN MATRIX INVERSION.
	INTEGER I,II,J,LL
	REAL SUM
	!                                   WHEN II IS SET TO A POSITIVE VALUE, IT WILL BECOME THE IN-
	II=0
	!                                         DEX OF THE ϬRST NONVANISHING ELEMENT OF B. WE NOW DO
	DO I=1,N
	!                                         THE FORWARD SUBSTITUTION, EQUATION (2.3.6). THE ONLY NEW
	LL=INDX(I)
	!                                         WRINKLE IS TO UNSCRAMBLE THE PERMUTATION AS WE GO.
	SUM=B(LL)
	B(LL)=B(I)
	IF (II.NE.0)THEN
	DO J=II,I-1
	SUM=SUM-A(I,J)*B(J)
	ENDDO
	ELSE IF (SUM.NE.0.) THEN
	!                                   A NONZERO ELEMENT WAS ENCOUNTERED, SO FROM NOW ON WE WILL
	II=I
	!                                         HAVE TO DO THE SUMS IN THE LOOP ABOVE.
	ENDIF
	B(I)=SUM
	ENDDO
	!                                   NOW WE DO THE BACKSUBSTITUTION, EQUATION (2.3.7).
	DO I=N,1,-1
	SUM=B(I)
	DO J=I+1,N
	SUM=SUM-A(I,J)*B(J)
	ENDDO
	!                                   STORE A COMPONENT OF THE SOLUTION VECTOR X.
	B(I)=SUM/A(I,I)
	ENDDO
	!                                   ALL DONE!
	RETURN
	END SUBROUTINE LUBKSB

      !.................................................................

	SUBROUTINE LUDCMP(A,N,DI,INDX,D)
        implicit none
	INTEGER N,DI,INDX(N),NMAX
	REAL D,A(DI,DI),TINY
	PARAMETER (NMAX=500,TINY=1.0E-20)! LARGEST EXPECTED N, AND A SMALL NUMBER.
	!     GIVEN A MATRIX A(1:N,1:N), WITH PHYSICAL DIMENSION NP BY NP, THIS ROUTINE REPLACES IT BY
	!     THE LU DECOMPOSITION OF A ROWWISE PERMUTATION OF ITSELF. A AND N ARE INPUT. A IS OUTPUT,
	!     ARRANGED AS IN EQUATION (2.3.14) ABOVE; INDX(1:N) IS AN OUTPUT VECTOR THAT RECORDS THE
	!     ROW PERMUTATION EϬECTED BY THE PARTIAL PIVOTING; D IS OUTPUT AS ±1 DEPENDING ON WHETHER
	!     THE NUMBER OF ROW INTERCHANGES WAS EVEN OR ODD, RESPECTIVELY. THIS ROUTINE IS USED IN
	!     COMBINATION WITH LUBKSB TO SOLVE LINEAR EQUATIONS OR INVERT A MATRIX.
	INTEGER I,IMAX,J,K
	!                                          VV STORES THE IMPLICIT SCALING OF EACH ROW.
	REAL AAMAX,DUM,SUM,VV(NMAX)
	!                                          NO ROW INTERCHANGES YET.
	D=1.
	!                                          LOOP OVER ROWS TO GET THE IMPLICIT SCALING INFORMA-
	DO I=1,N
	!                                               TION.
	AAMAX=0.
	DO J=1,N
	IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
	ENDDO
	IF (AAMAX.EQ.0.) STOP
        !SINGULAR MATRIX IN LUDCMP !NO NONZERO LARGEST ELEMENT.
	! SAVE THE SCALING.
	VV(I)=1./AAMAX
	ENDDO
	!                                          THIS IS THE LOOP OVER COLUMNS OF CROUTS METHOD.
	DO J=1,N
	!                                          THIS IS EQUATION (2.3.12) EXCEPT FOR I = J.
	DO I=1,J-1
	SUM=A(I,J)
	DO K=1,I-1
	SUM=SUM-A(I,K)*A(K,J)
	ENDDO
	A(I,J)=SUM
	ENDDO
	!                                          INITIALIZE FOR THE SEARCH FOR LARGEST PIVOT ELEMENT.
	AAMAX=0.
	!                                          THIS IS I = J OF EQUATION (2.3.12) AND I = J + 1 . . . N
	DO I=J,N
	!                                               OF EQUATION (2.3.13).
	SUM=A(I,J)
	DO K=1,J-1
	SUM=SUM-A(I,K)*A(K,J)
	ENDDO
	A(I,J)=SUM
	!                                          FIGURE OF MERIT FOR THE PIVOT.
	DUM=VV(I)*ABS(SUM)
	!                                          IS IT BETTER THAN THE BEST SO FAR?
	IF (DUM.GE.AAMAX) THEN
	IMAX=I
	AAMAX=DUM
	ENDIF
	ENDDO
	!                                          DO WE NEED TO INTERCHANGE ROWS?
	IF (J.NE.IMAX)THEN
	!                                          YES, DO SO...
	DO  K=1,N
	DUM=A(IMAX,K)
	A(IMAX,K)=A(J,K)
	A(J,K)=DUM
	ENDDO
	!                                          ...AND CHANGE THE PARITY OF D.
	D=-D
	!                                          ALSO INTERCHANGE THE SCALE FACTOR.
	VV(IMAX)=VV(J)
	ENDIF
	INDX(J)=IMAX
	IF(A(J,J).EQ.0.)A(J,J)=TINY
	!        IF THE PIVOT ELEMENT IS ZERO THE MATRIX IS SINGULAR (AT LEAST TO THE PRECISION OF THE AL-
	!        GORITHM). FOR SOME APPLICATIONS ON SINGULAR MATRICES, IT IS DESIRABLE TO SUBSTITUTE TINY
	!        FOR ZERO.
	!                                 NOW, ϬNALLY, DIVIDE BY THE PIVOT ELEMENT.
	IF(J.NE.N)THEN
	DUM=1./A(J,J)
	DO I=J+1,N
	A(I,J)=A(I,J)*DUM
	ENDDO
	ENDIF
	!                                 GO BACK FOR THE NEXT COLUMN IN THE REDUCTION.
	ENDDO
	RETURN
	END SUBROUTINE LUDCMP






end program kasuo_h
