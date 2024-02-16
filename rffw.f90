subroutine rffw(depth, vp, rho, vs, nlay, rayp, &
    gaussalp, delay, n, delta, rfcalc)
    implicit none
    integer, parameter :: maxdata=8192
    integer, intent(in) :: nlay
    integer,intent(in) :: n
    real(kind=4), intent(in), dimension(nlay) :: depth, vp, &
                                              rho, vs
    real(kind=4), intent(in) :: rayp, gaussalp, delay, delta
    real(kind=4), dimension(maxdata) :: rf
    integer :: rfdim
    real(kind=4), dimension((nlay*4)-1) :: modelin
    real(kind=4), dimension(n), intent(out) :: rfcalc
    real(kind=4), dimension(nlay) ::thick
    integer i, moddim

    do i=1, nlay-1
        if (i .eq. 1) then
            thick(i) =  depth(i)
            print*, thick(i), depth(i)
        else
            thick(i) = depth(i) - depth(i-1)
            print*, thick(i), depth(i), depth(i-1)
        end if
    enddo

    modelin(1:nlay) = vs
    modelin(nlay+1:nlay*2) = vp/vs
    modelin((nlay*2)+1:nlay*3) = rho
    modelin((nlay*3)+1:(nlay*4)-1) = thick(1:nlay-1)

    moddim = size(modelin)

!    do i = 1, (nlay*4)-1
!        write(*,*) i, modelin(i)
!    enddo
!    write(*,*) moddim

    call zero(rf, maxdata)

    call hrftn(modelin, moddim, rayp, &
       gaussalp, delay, n, delta, rf)

    rfdim = n
    rfcalc = rf(1:rfdim)

!     return
end subroutine rffw
