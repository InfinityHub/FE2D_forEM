 MODULE sparse_interface_mod
  IMPLICIT NONE

  ! Interface blocks for all Sparskit routines used by "g3dfd" and dependencies.

  INTERFACE

      ! From "formats.f" ...

      subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
        real*8 a(*),ao(*),x
        integer*4 nnz
      integer ir(*),jc(*),jao(*),iao(*)
      end

      subroutine csrcoo (nrow,job,nzmax,a,ja,ia,nnz,ao,ir,jc,ierr)
        real*8 a(*),ao(*) 
        integer ir(*),jc(*),ja(*),ia(nrow+1)
      end subroutine csrcoo

    
      subroutine csrcsc (n,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n+1),ja(*),jao(*)
      real*8  a(*),ao(*)
      end

      subroutine csrcsc2 (n,n2,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n2+1),ja(*),jao(*)
      real*8  a(*),ao(*)
      end

      ! From "blassm.f" ...

      subroutine amub (nrow,ncol,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr) 
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(*),ic(*),iw(ncol)
      integer*8  nzmax
      end

      subroutine aplb (nrow,ncol,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1),iw(ncol)
      end

      subroutine diamua (nrow,job, a, ja, ia, diag, b, jb, ib)
      real*8 a(*), b(*), diag(nrow), scal
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1) 
      end

      subroutine amudia (nrow,job, a, ja, ia, diag, b, jb, ib)
      real*8 a(*), b(*), diag(nrow) 
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1) 
      end

      subroutine apldia (nrow, job, a, ja, ia, diag, b, jb, ib, iw) 
      real*8 a(*), b(*), diag(nrow) 
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1), iw(*)
      end

      ! From "matvec.f" ...

      subroutine amux (n, x, y, a,ja,ia) 
      real*8  x(*), y(*), a(*) 
      integer n, ja(*), ia(*)
      end

      subroutine atmux (n, x, y, a, ja, ia)
      real*8 x(*), y(*), a(*) 
      integer n, ia(*), ja(*)
      end

      subroutine atmuxr (m, n, x, y, a, ja, ia)
      real*8 x(*), y(*), a(*) 
      integer m, n, ia(*), ja(*)
      end

      ! From "iters.f" ...

      subroutine bcgstab(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,8)
      end

      subroutine bcgstab_inv(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,8)
      end

      subroutine cg(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,*)
      end

      subroutine cgnr(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,*)
      end

      subroutine bcg(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,*)
      end

      !vvv For testing.
      subroutine dqgmres(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(*)
      end
      !^^^

      !vvv For testing.
      subroutine gmres(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(*)
      end
      !^^^

      !vvv For testing.
      subroutine gmres_inv(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(*)
      end
      !^^^

      ! From "ilut.f" ...

      subroutine ilutp(n,a,ja,ia,lfil,droptol,permtol,mbloc,alu,jlu,ju,iwk,w,jw,iperm,ierr)
      integer n,ja(*),ia(n+1),lfil,jlu(*),ju(n),jw(2*n),iwk,iperm(2*n),ierr
      real*8 a(*), alu(*), w(n), droptol
      integer k,i,j,jrow,ju0,ii,j1,j2,jpos,len,imax,lenu,lenl,jj,mbloc,icut
      real*8 s, tmp, tnorm,xmax,xmax0, fact, abs, t, permtol
      end

      subroutine ilut(n,a,ja,ia,lfil,droptol,alu,jlu,ju,iwk,w,jw,ierr)
      implicit none 
      integer n 
      real*8 a(*),alu(*),w(n),droptol
      integer ja(*),ia(n+1),jlu(*),ju(n),jw(2*n),lfil,iwk,ierr
      end

      subroutine lusol(n, y, x, alu, jlu, ju)
      real*8 x(n), y(n), alu(*)
      integer n, jlu(*), ju(*)
      end

      subroutine lutsol(n, y, x, alu, jlu, ju) 
      real*8 x(n), y(n), alu(*)
      integer n, jlu(*), ju(*)
      end

      ! From "unary.f" ...

      subroutine transp (nrow,ncol,a,ja,ia,iwk,ierr)
      integer nrow, ncol, ia(*), ja(*), iwk(*), ierr
      real*8 a(*) 
      end

      subroutine amubdg (nrow,ncol,ncolb,ja,ia,jb,ib,ndegr,nnz,iw) 
      integer ja(*),jb(*),ia(nrow+1),ib(ncol+1),ndegr(nrow),iw(ncolb) 
      end

      subroutine aplbdg (nrow,ncol,ja,ia,jb,ib,ndegr,nnz,iw) 
      integer ja(*),jb(*),ia(nrow+1),ib(nrow+1),iw(ncol),ndegr(nrow) 
      end


  END INTERFACE

 END MODULE sparse_interface_mod
