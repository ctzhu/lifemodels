!%%fortran
!-------------------------------------------------------------------------------
        subroutine gp3LLK(sx1, sx2, sy, sa, sb, sc, n)
        double precision, intent(in) :: sx1(n), sx2(n)
        double precision, intent(out) :: sy
        double precision :: sa, sb, sc
        real(kind=10) :: x1(n), x2(n)
        real(kind=10) :: y
        real(kind=10) :: ebx1, ebx2, a, b, c
        integer n, i

        x1 = sx1
        x2 = sx2
        a  = sa
        b  = sb
        c  = sc
        y  = 0.0_10
    
        do 100 i=1, n
            ebx1 = exp(b*x1(i))
            if (x1(i) == x2(i)) then
                y = y-(a/b)*(ebx1-1.0_10) &
                     -c*x1(i) &
                     +log(a*ebx1+c)
            else
                ebx2 = exp(b*x2(i))
                y = y+log(exp(-c*x1(i)-(a/b)*(ebx1-1.0_10)) &
                         -exp(-c*x2(i)-(a/b)*(ebx2-1.0_10)))
            end if
100     continue
        sy = dble(y)
        end subroutine gp3LLK

!Gompertz-Makeham distribution. 3 parameters
!Likelihood and neg-likelihood function with paramerters in log-scale 

        subroutine gp3LLK_lp(p, x, y, m, n)
        double precision, intent(in) :: x(2,n), p(m)        
        double precision, intent(out) :: y
        double precision :: x1(n), x2(n)
        double precision :: a, b, c
        integer :: n, m
        x1 = x(1, 1:n)
        x2 = x(2, 1:n)
        a = exp(p(1))
        b = exp(p(2))
        c = exp(p(3)) - a !TODO
        if ((p(3) < - 35.0d0) .or.  (p(3) > 5.0d0)) then
            call gp2LLK(x1, x2, y, a, b, n)
        else
            call gp3LLK(x1, x2, y, a, b, c, n)
        end if
        y = -y
        end subroutine gp3LLK_lp
        
!-------------------------------------------------------------------------------
        subroutine gl3LLK(sx1, sx2, sy, sa, sb, sd, n)
        double precision, intent(in) :: sx1(n), sx2(n)
        double precision, intent(out) :: sy
        double precision :: sa, sb, sd
        real(kind=10) :: x1(n), x2(n)
        real(kind=10) :: y
        real(kind=10) :: ebx1, ebx2, a, b, d, a_bd
        integer n, i

        x1   = sx1
        x2   = sx2
        a    = sa
        b    = sb
        d    = sd
        y    = 0.0_10
        a_bd = a/b/d
    
        do 100 i=1, n
            ebx1 = exp(b*x1(i))
            if (x1(i) == x2(i)) then
                y = y+a_bd*(log(d+1.0_10)-log(d*ebx1+1.0_10)) &
                     +log(a*ebx1)-log(1.0_10+d*ebx1)
            else
                ebx2 = exp(b*x2(i))
                y = y+log(exp(a_bd*(log(d+1.0_10)-log(d*ebx1+1.0_10))) &
                         -exp(a_bd*(log(d+1.0_10)-log(d*ebx2+1.0_10))))
        end if
        sy = dble(y)
100     continue
        end subroutine gl3LLK

!Gompertz-Logistic distribution. 3 parameters
!Likelihood and neg-likelihood function with paramerters in log-scale 

        subroutine gl3LLK_lp(p, x, y, m, n)
        double precision, intent(in) :: x(2,n), p(m)        
        double precision, intent(out) :: y
        double precision :: x1(n), x2(n)
        double precision :: a, b, d
        integer :: n, m
        x1 = x(1, 1:n)
        x2 = x(2, 1:n)
        a = exp(p(1))
        b = exp(p(2))
        d = exp(p(3))
        if ((p(3) < - 35.0d0) .or.  (p(3) > 5.0d0)) then
            call gp2LLK(x1, x2, y, a, b, n)
        else
            call gl3LLK(x1, x2, y, a, b, d, n)
        end if
        y = -y
        end subroutine gl3LLK_lp
        
!-------------------------------------------------------------------------------
        subroutine gl4LLK(sx1, sx2, sy, sa, sb, sc, sd, n)
        double precision, intent(in)  :: sx1(n), sx2(n)
        double precision, intent(out) :: sy
        double precision :: sa, sb, sc, sd
        integer :: n, i

        real(kind=10) :: x1(n), x2(n)
        real(kind=10) :: y
        real(kind=10) :: ebx1, ebx2, a, b, c, d, a_bd

        x1 = sx1
        x2 = sx2
        a  = sa
        b  = sb
        c  = sc
        d  = sd
        y = 0.0_10
        a_bd = a/b/d
    
        do 100 i=1, n
            ebx1 = exp(b*x1(i))
            if (x1(i) == x2(i)) then
                y = y+a_bd*(log(d+1.0_10)-log(d*ebx1+1.0_10)) &
                     -c*x1(i) &
                     +log(a*ebx1/(1.0_10+d*ebx1)+c)
            else
                ebx2 = exp(b*x2(i))
                y = y+log(exp(a_bd*(log(d+1.0_10)-log(d*ebx1+1.0_10)) &
                               -c*x1(i)) &
                         -exp(a_bd*(log(d+1.0_10)-log(d*ebx2+1.0_10)) &
                               -c*x2(i)))
            end if
        sy = dble(y)
100     continue
        end subroutine gl4LLK

!Gompertz-Logistic-Makeham distribution. 4 parameters
!Likelihood and neg-likelihood function with paramerters in log-scale 

        subroutine gl4LLK_lp(p, x, y, m, n)
        double precision, intent(in) :: x(2,n), p(m)        
        double precision, intent(out) :: y
        double precision :: x1(n), x2(n)
        double precision :: a, b, c, d
        integer :: n, m
        x1 = x(1, 1:n)
        x2 = x(2, 1:n)
        a = exp(p(1))
        b = exp(p(2))
        d = exp(p(4))
        c = exp(p(3)) - a/(1.0d0+d) !TODO
        if      (((p(3) < -50.0d0) .or.  (p(3) > 5.0d0)) &
          .and.  ((p(4) > -50.0d0) .and. ((p(4)-p(1)) < 5.0d0)))then
            call gl3LLK(x1, x2, y, a, b, d, n)
        else if (((p(3) > -50.0d0) .and. (p(3) < 5.0d0)) &
          .and.  ((p(4) < -50.0d0) .or.  ((p(4)-p(1)) > 5.0d0)))then
            call gp3LLK(x1, x2, y, a, b, c, n)
        else if (((p(3) < -50.0d0) .or.  (p(3) > 5.0d0)) &
          .and.  ((p(4) < -50.0d0) .or.  ((p(4)-p(1)) > 5.0d0)))then
            call gp2LLK(x1, x2, y, a, b, n)
        else
            call gl4LLK(x1, x2, y, a, b, c, d, n)
        end if
        y = -y
        end subroutine gl4LLK_lp

!-------------------------------------------------------------------------------
        subroutine gp2LLK(x1, x2, y, a, b, n)
        double precision, intent(in) :: x1(n), x2(n)
        double precision, intent(out) :: y
        double precision :: ebx1, ebx2, a, b
        integer n, i
        y = 0.0d+0
    
        do 100 i=1, n
            ebx1 = exp(b*x1(i))
            if (x1(i) == x2(i)) then
                y = y-a/b*(ebx1-1.0d0) &
                     +log(a*ebx1)
            else
                ebx2 = exp(b*x2(i))
                y = y+log(exp(-a/b*(ebx1-1.0d0)) &
                         -exp(-a/b*(ebx2-1.0d0)))
            end if
100     continue
        end subroutine gp2LLK

!Gompertz distribution. 2 parameters
!Likelihood and neg-likelihood function with paramerters in log-scale 

        subroutine gp2LLK_lp(p, x, y, m, n)
        double precision, intent(in) :: x(2,n), p(m)        
        double precision, intent(out) :: y
        double precision :: x1(n), x2(n)
        double precision :: a, b
        integer :: n, m
        x1 = x(1, 1:n)
        x2 = x(2, 1:n)
        a = exp(p(1))
        b = exp(p(2))
        call gp2LLK(x1, x2, y, a, b, n)
        y = -y
        end subroutine gp2LLK_lp

!-------------------------------------------------------------------------------
        subroutine wb2LLK(x1, x2, y, l, k, n)
        double precision, intent(in) :: x1(n), x2(n)
        double precision, intent(out) :: y
        double precision :: l, k
        integer n, i
        y = 0.0d+0
    
        do 100 i=1, n
            if (x1(i) == x2(i)) then
                if (x1(i)==0.0d0) then
                    y = y-0.0d0 !pdf is 0 at this point, what TODO?
                else
                    y = y-(x1(i)/l)**k &
                         +log(k/l)+(k-1.0d0)*log(x1(i)/l)
                end if
            else
                y = y+log(exp(-(x1(i)/l)**k)          &
                         -exp(-(x2(i)/l)**k))
            end if
100     continue
        end subroutine wb2LLK

!Weibull distribution. 2 parameters
!Likelihood and neg-likelihood function with paramerters in log-scale 

        subroutine wb2LLK_lp(p, x, y, m, n)
        double precision, intent(in) :: x(2,n), p(m)        
        double precision, intent(out) :: y
        double precision :: x1(n), x2(n)
        double precision :: l, k
        integer :: n, m
        x1 = x(1, 1:n)
        x2 = x(2, 1:n)
        l = exp(p(1))
        k = exp(p(2))
        call wb2LLK(x1, x2, y, l, k, n)
        y = -y
        end subroutine wb2LLK_lp

!-------------------------------------------------------------------------------
        subroutine wb3LLK(x1, x2, y, l, k, c, n)
        double precision, intent(in) :: x1(n), x2(n)
        double precision, intent(out) :: y
        double precision :: l, k, c
        integer n, i
        y = 0.0d+0
    
        do 100 i=1, n
            if (x1(i) == x2(i)) then
                if (x1(i) == 0.0d0) then
                    y = y-log(c) !pdf=c, as sf=1 and hzd=c at 0
                else
                    y = y-c*x1(i)-(x1(i)/l)**k &
                         +log((k/l)*((x1(i)/l)**(k-1.0d0))+c)
                end if
            else
                y = y+log(exp(-c*x1(i)-(x1(i)/l)**k) &
                         -exp(-c*x2(i)-(x2(i)/l)**k))
            end if
100     continue
        end subroutine wb3LLK

!Weibull-Makeham distribution. 3 parameters
!Likelihood and neg-likelihood function with paramerters in log-scale 

        subroutine wb3LLK_lp(p, x, y, m, n)
        double precision, intent(in) :: x(2,n), p(m)        
        double precision, intent(out) :: y
        double precision :: x1(n), x2(n)
        double precision :: l, k, c
        integer :: n, m
        x1 = x(1, 1:n)
        x2 = x(2, 1:n)
        l = exp(p(1))
        k = exp(p(2))
        c = exp(p(3))
        if ((p(3) < - 35.0d0) .or.  (p(3) > 5.0d0)) then
            call wb2LLK(x1, x2, y, l, k, n)
        else
            call wb3LLK(x1, x2, y, l, k, c, n)
        end if
        y = -y
        end subroutine wb3LLK_lp

!-------------------------------------------------------------------------------
        subroutine lg2LLK(x1, x2, y, u, s, n)
        double precision, intent(in) :: x1(n), x2(n)
        double precision, intent(out) :: y
        double precision :: v1, v2, u, s
        integer n, i
        y = 0.0d+0
            
        do 100 i=1, n
            v1 = exp((u-x1(i))/s)
            if (x1(i) == x2(i)) then
                y = y+log(v1)-log(s)-2.0d0*log(1.0d0+v1)
            else
                v2 = exp((u-x2(i))/s)
                y = y+log(1.0d0/(1.0d0+v2) &
                         -1.0d0/(1.0d0+v1))
            end if
100     continue
        end subroutine lg2LLK

!Logistic distribution. 2 parameters
!Likelihood and neg-likelihood function with paramerters in log-scale 

        subroutine lg2LLK_lp(p, x, y, m, n)
        double precision, intent(in) :: x(2,n), p(m)        
        double precision, intent(out) :: y
        double precision :: x1(n), x2(n)
        double precision :: u, s
        integer :: n, m
        x1 = x(1, 1:n)
        x2 = x(2, 1:n)
        u = exp(p(1))
        s = exp(p(2))
        call lg2LLK(x1, x2, y, u, s, n)
        y = -y
        end subroutine lg2LLK_lp

!-------------------------------------------------------------------------------
        subroutine lg3LLK(sx1, sx2, sy, su, ss, sc, n)
        double precision, intent(in) :: sx1(n), sx2(n)
        double precision, intent(out) :: sy
        double precision :: su, ss, sc
        integer n, i

        real(kind=10) :: x1(n), x2(n)
        real(kind=10) :: y
        real(kind=10) :: v1, v2, u, s, c

        x1= sx1
        x2= sx2
        u = su
        s = ss
        c = sc
        y = 0.0_10
    
        do 100 i=1, n
            v1 = exp((u-x1(i))/s)
            if (x1(i) == x2(i)) then
                y = y+log(v1)-log(1.0_10+v1) &
                     -c*x1(i) &
                     +log(1.0_10/s/(1.0_10+v1)+c)
            else
                v2 = exp((u-x2(i))/s)
                y = y+log(v1/(1.0_10+v1)/exp(c*x1(i)) &
                         -v2/(1.0_10+v2)/exp(c*x2(i)))
            end if
100     continue
        sy = dble(y)
        end subroutine lg3LLK

!logistic-Makeham distribution. 3 parameters
!Likelihood and neg-likelihood function with paramerters in log-scale 

        subroutine lg3LLK_lp(p, x, y, m, n)
        double precision, intent(in) :: x(2,n), p(m)        
        double precision, intent(out) :: y
        double precision :: x1(n), x2(n)
        double precision :: u, s, c
        integer :: n, m
        x1 = x(1, 1:n)
        x2 = x(2, 1:n)
        u = exp(p(1))
        s = exp(p(2))
        c = exp(p(3)) - 1.0d0/s/(1.0d0+exp(u/s)) !TODO
        if ((p(3) < - 35.0d0) .or.  (p(3) > 5.0d0)) then
            call lg2LLK(x1, x2, y, u, s, n)
        else
            call lg3LLK(x1, x2, y, u, s, c, n)
        end if
        y = -y
        end subroutine lg3LLK_lp

!-------------------------------------------------------------------------------
        subroutine llg2LLK(x1, x2, y, a, b, n)
        double precision, intent(in) :: x1(n), x2(n)
        double precision, intent(out) :: y
        double precision :: v1, v2, a, b
        integer n, i
        y = 0.0d+0
    
        do 100 i=1, n
            v1 = (x1(i)/a)**b
            if (x1(i) == x2(i)) then
                if (x1(i)==0.0d0) then
                    y = y+0.0d0 !pdf is 0 at this point, what to do?
                else
                    y = y+log(b/x1(i))+log(v1)-2.0d0*log(1.0d0+v1)
                end if
            else
                v2 = (x2(i)/a)**b
                y = y+log(1.0d0/(1.0d0+v1)      &
                         -1.0d0/(1.0d0+v2))
            end if
100     continue
        end subroutine llg2LLK

!Log-logistic distribution. 2 parameters
!Likelihood and neg-likelihood function with paramerters in log-scale 

        subroutine llg2LLK_lp(p, x, y, m, n)
        double precision, intent(in) :: x(2,n), p(m)        
        double precision, intent(out) :: y
        double precision :: x1(n), x2(n)
        double precision :: a, b
        integer :: n, m
        x1 = x(1, 1:n)
        x2 = x(2, 1:n)
        a = exp(p(1))
        b = exp(p(2))
        call llg2LLK(x1, x2, y, a, b, n)
        y = -y
        end subroutine llg2LLK_lp

!-------------------------------------------------------------------------------
        subroutine llg3LLK(x1, x2, y, a, b, c, n)
        double precision, intent(in) :: x1(n), x2(n)
        double precision, intent(out) :: y
        double precision :: v1, v2, a, b, c, b_a
        integer n, i
        y = 0.0d+0
        b_a = b/a
    
        do 100 i=1, n
            v1 = (x1(i)/a)**(-b)
            if (x1(i) == x2(i)) then
                if (x1(i) == 0.0d0) then
                    y=y-log(c) !pdf=c, as sf=1 and hzd=c here.
                else
                    y=y+log(v1)-log(1.0d0+v1) &
                       -c*x1(i) &
                       +log(b/x1(i)/(1.0d0+v1)+c)
                end if
            else
                v2 = (x2(i)/a)**(-b)
                y = y+log(v1/(1.0d0+v1)/exp(x1(i)*c)   &
                         -v2/(1.0d0+v2)/exp(x2(i)*c))
            end if
100     continue
        end subroutine llg3LLK

!Log-logistic-Makeham distribution. 3 parameters
!Likelihood and neg-likelihood function with paramerters in log-scale 

        subroutine llg3LLK_lp(p, x, y, m, n)
        double precision, intent(in) :: x(2,n), p(m)        
        double precision, intent(out) :: y
        double precision :: x1(n), x2(n)
        double precision :: a, b, c
        integer :: n, m
        x1 = x(1, 1:n)
        x2 = x(2, 1:n)
        a = exp(p(1))
        b = exp(p(2))
        c = exp(p(3))
        if ((p(3) < - 35.0d0) .or.  (p(3) > 5.0d0)) then
            call llg2LLK(x1, x2, y, a, b, n)
        else
            call llg3LLK(x1, x2, y, a, b, c, n)
        end if
        y = -y
        end subroutine llg3LLK_lp

!-------------------------------------------------------------------------------
        subroutine scan_m(x, ly, lp, func, k, n, m)
        double precision, intent(in) :: x(2, n), lp(k, m)
        double precision, intent(out):: ly(k)
        integer :: n, m, k, i
        interface afunc
            subroutine func(p, x, y, m, n)
                double precision, dimension(2,n) :: x
                double precision, dimension(m) :: p
                double precision, intent (out) :: y
                integer :: m, n
            end subroutine func
        end interface afunc
        do i=1, k
            call func(lp(i, 1:m), x, ly(i), m, n)
            !Ideally: func should be ??LLk_lp() subroutines
        end do
        end subroutine scan_m

        subroutine scan_2(x, yL, aL, bL, n, m)
        double precision, intent(in) :: x(2,n), aL(m), bL(m)
        double precision, intent(out) :: yL(m)
        double precision :: x1(n), x2(n)
        integer :: n, m, i
        x1 = x(1, 1:n)
        x2 = x(2, 1:n)        
        do 101 i=1, m
            call gp2LLK(x1, x2, yL(i), aL(i), bL(i), n)
101     continue
        end subroutine scan_2
    
        subroutine scan_range_2(x, ar, br, m, yL, aL, bL, n)
        double precision, intent(in) :: x(2,n), ar(2), br(2)
        double precision, intent(out) :: yL(m*m), aL(m*m), bL(m*m)
        double precision :: x1(n), x2(n)
        double precision :: a_low, b_low, a_step, b_step
        integer :: n, m, i, j, k
        x1 = x(1, 1:n)
        x2 = x(2, 1:n)
        a_low = ar(1)
        b_low = br(1)
        a_step = (ar(2)-ar(1))/m
        b_step = (br(2)-br(1))/m
        
        do 102 i=0, m-1
            do j=1, m
                k = i*m + j
                aL(k) = a_low + i*a_step
                bL(k) = b_low + j*b_step
                call gp2LLK(x1, x2, yL(k), exp(aL(k)), exp(bL(k)), n)
            end do
102     continue
        end subroutine scan_range_2

!-------------------------------------------------------------------------------
!Using REAL(10) in internal.
!%%fortran
     subroutine derivative(p, x, r, func, m, n)
       double precision, intent (in) :: x(2,n), p(m)
       double precision, intent(out) :: r(m)
       double precision :: dp(m)
       double precision :: d, y ,dy, scl_fct
       integer :: m, n
       interface afunc
          subroutine func(p, x, y, m, n)
             double precision, dimension(2,n) :: x
             double precision, dimension(m) :: p
             double precision, intent (out) :: y
             integer :: m, n
          end subroutine func
       end interface afunc
       d = 0.000001d0
       call func(p, x, y, m, n)
       if (abs(y)>100000.0d0) then
          scl_fct = 0.0000001d0
       else if (abs(y)<0.00001d0) then
          scl_fct = 10000000.0d0
       else
          scl_fct = 1.0d0
       end if
       do i=1, m
          dp = p
          dp(i) = dp(i) + d
          call func(dp, x, dy, m, n)
          r(i) = dy*scl_fct-y*scl_fct
       end do
       r = r/d/scl_fct
       end subroutine derivative

     subroutine derivative2(p, x, r, func, m, n)
       double precision, intent (in) :: x(2,n), p(m)
       double precision, intent(out) :: r(m, m)
       double precision :: dp(m), dy(m), y(m)
       double precision :: d
       integer :: m, n
       interface afunc
          subroutine func(p, x, y, m, n)
             double precision, dimension(2,n) :: x
             double precision, dimension(m) :: p
             double precision, intent (out) :: y
             integer :: m, n
          end subroutine func
       end interface afunc
       d = 0.000001d0
       call derivative(p, x, y, func, m, n)
       do i=1, m
          dp = p
          dp(i) = dp(i) + d
          call derivative(dp, x, dy, func, m, n)
          r(i, 1:m) = (dy-y)/d
       end do
       end subroutine derivative2

     subroutine ederivative(p, x, r, func, m, n)
       double precision, intent (in) :: x(2,n), p(m)
       double precision, intent(out) :: r(m)
       double precision :: dp(m)
       double precision :: d, y ,dy
       real(10) :: pd, py, pdy, scl_fct !to map to C long double
       integer  :: m, n
       interface afunc
          subroutine func(p, x, y, m, n)
             double precision, dimension(2,n) :: x
             double precision, dimension(m) :: p
             double precision, intent (out) :: y
             integer :: m, n
          end subroutine func
       end interface afunc
       d  = 0.000001d0
       pd = 0.000001_10
       call func(p, x, y, m, n)
       if (abs(y)>100000.0d0) then
          scl_fct = 0.0000001_10
       else if (abs(y)<0.00001d0) then
          scl_fct = 10000000.0_10
       else
          scl_fct = 1.0_10
       end if
       py = 0.0_10 + y
       do i=1, m
          dp = p
          dp(i) = dp(i) + d
          call func(dp, x, dy, m, n)
          pdy  = 0.0_10 + dy
          r(i) = dble((pdy*scl_fct-py*scl_fct)/scl_fct) !The only step in REAL(10)
          !r(i) = (dy-y)/d
       end do
       r = r/d
       end subroutine ederivative

     subroutine ederivative2(p, x, r, func, m, n)
       double precision, intent (in) :: x(2,n), p(m)
       double precision, intent(out) :: r(m, m)
       double precision :: dp(m), dy(m), y(m)
       double precision :: d
       real(10) :: pd, py(m), pdy(m) 
       integer  :: m, n
       interface afunc
          subroutine func(p, x, y, m, n)
             double precision, dimension(2,n) :: x
             double precision, dimension(m) :: p
             double precision, intent (out) :: y
             integer :: m, n
          end subroutine func
       end interface afunc
       d  = 0.000001d0
       pd = 0.000001_10
       call derivative(p, x, y, func, m, n)
       py = 0.0_10 + y 
       do i=1, m
          dp = p
          dp(i) = dp(i) + d
          call derivative(dp, x, dy, func, m, n)
          pdy = 0.0_10 + dy
          r(i, 1:m) = dble((pdy-py)/pd) !The only step in REAL(10)
          !r(i, 1:m) = (dy-y)/d
       end do
       end subroutine ederivative2

!do not need to change (2,n) in subrt derivative
!to (n) when calling derviative on this trial function
!as x has (n) dimemsion here, the extra dimension in (2,n)
!is not going to do anything, anyway. 
!    subroutine unifunc(p, x, y, m, n)
!       double precision, intent (in) :: x(n), p(m)
!       double precision, intent(out) :: y
!       integer :: m, n
!       y = x(1)*p(1)*p(1)+x(2)*p(2)
!       end subroutine unifunc
                                
!    subroutine f(x, r, m)
!        double precision, intent(out):: r(m)
!        double precision, intent (in) :: x(2,m)
!        integer :: m
!        r = x(1, 1:m)
!    end subroutine f


!-------------------------------------------------------------------------------
      subroutine inverse(m,c,n)
      !============================================================
      ! Inverse matrix
      ! Method: Based on Doolittle LU factorization for Ax=b
      ! Alex G. December 2009
      !-----------------------------------------------------------
      ! input ...
      ! a(n,n) - array of coefficients for matrix A
      ! n      - dimension
      ! output ...
      ! c(n,n) - inverse matrix of A
      ! comments ...
      ! the original matrix a(n,n) will be destroyed 
      ! during the calculation
      !-----------------------------------------------------------
      ! Modified for numpy/f2py use.
      ! Matrix 'a' is created from 'm' to avoid altering the input
      !   array, a precaution for using with f2py.
      ! C-T Zhu, December 2014
      !===========================================================
      implicit none 
      double precision, intent(in)  :: m(n,n)
      double precision, intent(out) :: c(n,n)
      double precision :: a(n,n), L(n,n), U(n,n), b(n), d(n), x(n)
      double precision              :: coeff
      integer                       :: i, j, k, n
      
      ! step 0: initialization for matrices L and U and b
      ! Fortran 90/95 aloows such operations on matrices
      L = 0.0d0
      U = 0.0d0
      b = 0.0d0
      a = m
      
      ! step 1: forward elimination
      do k=1, n-1
          do i=k+1,n
              coeff  = a(i,k)/a(k,k)
              L(i,k) = coeff
              do j=k+1,n
                  a(i,j) = a(i,j)-coeff*a(k,j)
              end do
          end do
      end do
      
      ! Step 2: prepare L and U matrices 
      ! L matrix is a matrix of the elimination coefficient
      ! + the diagonal elements are 1.0
      do i=1,n
          L(i,i) = 1.0d0
      end do
      ! U matrix is the upper triangular part of A
      do j=1,n
          do i=1,j
              U(i,j) = a(i,j)
          end do
      end do
      
      ! Step 3: compute columns of the inverse matrix C
      do k=1,n
          b(k)=1.0d0
          d(1) = b(1)
        ! Step 3a: Solve Ld=b using the forward substitution
          do i=2,n
              d(i) = b(i)
              do j=1,i-1
                  d(i) = d(i) - L(i,j)*d(j)
              end do
          end do
        ! Step 3b: Solve Ux=d using the back substitution
          x(n)=d(n)/U(n,n)
          do i = n-1,1,-1
              x(i) = d(i)
              do j=n,i+1,-1
                  x(i) = x(i)-U(i,j)*x(j)
              end do
              x(i) = x(i)/u(i,i)
          end do
        ! Step 3c: fill the solutions x(n) into column k of C
          do i=1,n
              c(i,k) = x(i)
          end do
          b(k) = 0.0d0
      end do
      end subroutine inverse


!-------------------------------------------------------------------------------
      subroutine einverse(m,r,n)
      !============================================================
      ! Inverse matrix
      ! Method: Based on Doolittle LU factorization for Ax=b
      ! Alex G. December 2009
      !-----------------------------------------------------------
      ! input ...
      ! a(n,n) - array of coefficients for matrix A
      ! n      - dimension
      ! output ...
      ! c(n,n) - inverse matrix of A
      ! comments ...
      ! the original matrix a(n,n) will be destroyed 
      ! during the calculation
      !-----------------------------------------------------------
      ! Modified for numpy/f2py use, long_float internals.
      ! Matrix 'a' is created from 'm' to avoid altering the input
      !   array, a precaution for using with f2py.
      ! C-T Zhu, December 2014
      !===========================================================
      implicit none 
      double precision, intent(in)  :: m(n,n)
      double precision, intent(out) :: r(n,n)
      real(kind=10)  :: a(n,n), c(n,n), L(n,n), U(n,n), b(n), d(n), x(n)
      real(kind=10)                 :: coeff
      integer                       :: i, j, k, n
      
      ! step 0: initialization for matrices L and U and b
      ! Fortran 90/95 aloows such operations on matrices
      L = 0.0_10
      U = 0.0_10
      b = 0.0_10
      a = m+0.0_10
      
      ! step 1: forward elimination
      do k=1, n-1
          do i=k+1,n
              coeff  = a(i,k)/a(k,k)
              L(i,k) = coeff
              do j=k+1,n
                  a(i,j) = a(i,j)-coeff*a(k,j)
              end do
          end do
      end do
      
      ! Step 2: prepare L and U matrices 
      ! L matrix is a matrix of the elimination coefficient
      ! + the diagonal elements are 1.0
      do i=1,n
          L(i,i) = 1.0d0
      end do
      ! U matrix is the upper triangular part of A
      do j=1,n
          do i=1,j
              U(i,j) = a(i,j)
          end do
      end do
      
      ! Step 3: compute columns of the inverse matrix C
      do k=1,n
          b(k)=1.0d0
          d(1) = b(1)
        ! Step 3a: Solve Ld=b using the forward substitution
          do i=2,n
              d(i)=b(i)
              do j=1,i-1
                  d(i) = d(i) - L(i,j)*d(j)
              end do
          end do
        ! Step 3b: Solve Ux=d using the back substitution
          x(n)=d(n)/U(n,n)
          do i = n-1,1,-1
              x(i) = d(i)
              do j=n,i+1,-1
                  x(i) = x(i)-U(i,j)*x(j)
              end do
              x(i) = x(i)/u(i,i)
          end do
        ! Step 3c: fill the solutions x(n) into column k of C
          do i=1,n
              c(i,k) = x(i)
          end do
          b(k) = 0.0d0
      end do
      !print 
      r = dble(c)
      end subroutine einverse


