!-------------------------------------------------------------------------------
        subroutine gp2_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: a, b
        double precision :: ebx(n)
        integer :: m, n
        a = p(1)
        b = p(2)
        ebx = exp(b*x)
        hzd = a*ebx
        sf  = exp(-a/b*(ebx-1.0d0))
        pdf = hzd*sf
        cdf = 1.0d0-sf
        end subroutine gp2_cpsh

!Gompertz distribution. 2 parameters
!Return cdf, pdf, sf and hazard function, given parameter P and data X
!The 'name'lg_cpsh() takes parameter in log-scale.

        subroutine gp2lg_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: expp(m)
        integer :: n, m
        expp = exp(p)
        call gp2_cpsh(expp, x, cdf, pdf, sf, hzd, m, n)
        end subroutine gp2lg_cpsh

!-------------------------------------------------------------------------------
        subroutine gp3_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: a, b, c
        double precision :: ebx(n)
        integer :: m, n
        a = p(1)
        b = p(2)
        c = p(3)
        ebx = exp(b*x)
        hzd = a*ebx+c
        sf  = exp(-c*x-a/b*(ebx-1.0d0))
        pdf = hzd*sf
        cdf = 1.0d0-sf
        end subroutine gp3_cpsh

!Gompertz-Makeham distribution. 3 parameters
!Return cdf, pdf, sf and hazard function, given parameter P and data X
!The 'name'lg_cpsh() takes parameter in log-scale.

        subroutine gp3lg_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: expp(m)
        integer :: n, m
        expp = exp(p)
        expp(3) = expp(3)-expp(1) !TODO
        call gp3_cpsh(expp, x, cdf, pdf, sf, hzd, m, n)
        end subroutine gp3lg_cpsh

!-------------------------------------------------------------------------------
        subroutine gl3_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: a, b, d
        double precision :: ebx(n)
        integer :: m, n
        a = p(1)
        b = p(2)
        d = p(3)
        ebx = exp(b*x)
        hzd = a*ebx/(1.0d0+d*ebx)
        sf  = exp(a*log((d+1.0d0)/(d*ebx+1.0d0))/b/d)
        pdf = hzd*sf
        cdf = 1.0d0-sf
        end subroutine gl3_cpsh

!Gompertz-logistic distribution. 3 parameters
!Return cdf, pdf, sf and hazard function, given parameter P and data X
!The 'name'lg_cpsh() takes parameter in log-scale.

        subroutine gl3lg_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: expp(m)
        integer :: n, m
        expp = exp(p)
        call gl3_cpsh(expp, x, cdf, pdf, sf, hzd, m, n)
        end subroutine gl3lg_cpsh

!-------------------------------------------------------------------------------
        subroutine gl4_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: a, b, c, d
        double precision :: ebx(n)
        integer :: m, n
        a = p(1)
        b = p(2)
        c = p(3) 
        d = p(4)
        ebx = exp(b*x)
        hzd = a*ebx/(1.0d0+d*ebx)+c
        sf  = exp(a*log((d+1.0d0)/(d*ebx+1.0d0))/b/d-c*x)
        pdf = hzd*sf
        cdf = 1.0d0-sf
        end subroutine gl4_cpsh

!Gompertz-logistic-Makeham distribution. 4 parameters
!Return cdf, pdf, sf and hazard function, given parameter P and data X
!The 'name'lg_cpsh() takes parameter in log-scale.

        subroutine gl4lg_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: expp(m)
        integer :: n, m
        expp = exp(p)
        expp(3) = expp(3) - expp(1)/(1.0d0+expp(4)) !TODO
        call gl4_cpsh(expp, x, cdf, pdf, sf, hzd, m, n)
        end subroutine gl4lg_cpsh

!-------------------------------------------------------------------------------
        subroutine wb2_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: l, k
        double precision :: v(n)
        integer :: m, n
        l = p(1)
        k = p(2)
        v = (x/l)**(k-1.0d0) !avoid devide by x (zero-division error)
        hzd = v*k/l
        sf  = exp(-v*x/l)
        pdf = hzd*sf
        cdf = 1.0d0-sf
        end subroutine wb2_cpsh

!Weibull distribution. 2 parameters
!Return cdf, pdf, sf and hazard function, given parameter P and data X
!The 'name'lg_cpsh() takes parameter in log-scale.

        subroutine wb2lg_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: expp(m)
        integer :: n, m
        expp = exp(p)
        call wb2_cpsh(expp, x, cdf, pdf, sf, hzd, m, n)
        end subroutine wb2lg_cpsh

!-------------------------------------------------------------------------------
        subroutine wb3_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: c, l, k
        double precision :: v(n)
        integer :: m, n
        c = p(3)
        l = p(1)
        k = p(2)
        v = (x/l)**(k-1.0d0)
        hzd = v*k/l+c
        sf  = exp(-v*x/l-c*x)
        pdf = hzd*sf
        cdf = 1.0d0-sf
        end subroutine wb3_cpsh

!Weibull-Makeham distribution. 3 parameters
!Return cdf, pdf, sf and hazard function, given parameter P and data X
!The 'name'lg_cpsh() takes parameter in log-scale.

        subroutine wb3lg_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: expp(m)
        integer :: n, m
        expp = exp(p)
        call wb3_cpsh(expp, x, cdf, pdf, sf, hzd, m, n)
        end subroutine wb3lg_cpsh

!-------------------------------------------------------------------------------
        subroutine lg2_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: s, u
        double precision :: v(n)
        integer :: m, n
        u = p(1)
        s = p(2)
        v = exp((u-x)/s)
        hzd = 1.0d0/s/(1.0d0+v)
        sf  = v/(1.0d0+v)
        pdf = hzd*sf
        cdf = 1.0d0-sf
        end subroutine lg2_cpsh

!Logistic distribution. 2 parameters
!Return cdf, pdf, sf and hazard function, given parameter P and data X
!The 'name'lg_cpsh() takes parameter in log-scale.

        subroutine lg2lg_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: expp(m)
        integer :: n, m
        expp = exp(p)
        call lg2_cpsh(expp, x, cdf, pdf, sf, hzd, m, n)
        end subroutine lg2lg_cpsh

!-------------------------------------------------------------------------------
        subroutine lg3_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: c, s, u
        double precision :: v(n)
        integer :: m, n
        c = p(3) 
        u = p(1)
        s = p(2)
        v = exp((u-x)/s)
        hzd = 1.0d0/s/(1.0d0+v)+c
        sf  = v/(1.0d0+v)/exp(c*x)
        pdf = hzd*sf
        cdf = 1.0d0-sf
        end subroutine lg3_cpsh

!Logistic-Makeham distribution. 3 parameters
!Return cdf, pdf, sf and hazard function, given parameter P and data X
!The 'name'lg_cpsh() takes parameter in log-scale.

        subroutine lg3lg_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: expp(m)
        integer :: n, m
        expp = exp(p)
        expp(3) = expp(3) - 1.0d0/expp(2)/(1.0d0+exp(expp(1)/expp(2))) !TODO        
        call lg3_cpsh(expp, x, cdf, pdf, sf, hzd, m, n)
        end subroutine lg3lg_cpsh

!-------------------------------------------------------------------------------
        subroutine llg2_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: a, b
        double precision :: v(n)
        integer :: m, n
        a = p(1)
        b = p(2)
        v = (x/a)**(b-1)
        hzd = b*v/a/(1.0d0+v*x/a)
        sf  = 1.0d0/(1.0d0+v*x/a)
        pdf = hzd*sf
        cdf = 1.0d0-sf
        end subroutine llg2_cpsh

!Log-logistic distribution. 2 parameters
!Return cdf, pdf, sf and hazard function, given parameter P and data X
!The 'name'lg_cpsh() takes parameter in log-scale.

        subroutine llg2lg_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: expp(m)
        integer :: n, m
        expp = exp(p)
        call llg2_cpsh(expp, x, cdf, pdf, sf, hzd, m, n)
        end subroutine llg2lg_cpsh

!-------------------------------------------------------------------------------
        subroutine llg3_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: a, b, c
        double precision :: v(n)
        integer :: m, n
        a = p(1)
        b = p(2)
        c = p(3)
        v = (x/a)**(b-1)
        hzd = b*v/a/(1.0d0+v*x/a)+c
        sf  = 1.0d0/(1.0d0+v*x/a)/exp(c*x)
        pdf = hzd*sf
        cdf = 1.0d0-sf
        end subroutine llg3_cpsh

!Log-logistic-Makehame distribution. X parameters
!Return cdf, pdf, sf and hazard function, given parameter P and data X
!The 'name'lg_cpsh() takes parameter in log-scale.

        subroutine llg3lg_cpsh(p, x, cdf, pdf, sf, hzd, m, n)
        double precision, intent(in) :: p(m), x(n)
        double precision, intent(out):: cdf(n), pdf(n), sf(n), hzd(n)
        double precision :: expp(m)
        integer :: n, m
        expp = exp(p)
        call llg3_cpsh(expp, x, cdf, pdf, sf, hzd, m, n)
        end subroutine llg3lg_cpsh
