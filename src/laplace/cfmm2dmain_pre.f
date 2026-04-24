cc Copyright (C) 2018-2019: Leslie Greengard, Zydrunas Gimbutas,
cc and Manas Rachh
cc Contact: greengard@cims.nyu.edu
cc
cc convert to Cauchy FMM - Travis Askham 2021/07/07
cc
cc Pre-computed list variant - plan-based interface 2026
cc  cfmm2dmain_pre: identical to cfmm2dmain but accepts pre-computed
cc  carray and interaction lists as arguments rather than computing them
cc  internally. Used by cfmm2d_execute_plan.
cc
cc This program is free software; you can redistribute it and/or modify
cc it under the terms of the GNU General Public License as published by
cc the Free Software Foundation; either version 2 of the License, or
cc (at your option) any later version.

      subroutine cfmm2dmain_pre(nd,eps,
     $     nsource,sourcesort,
     $     ifcharge,chargesort,
     $     ifdipole,dipstrsort,
     $     ntarget,targetsort,nexpc,expcsort,
     $     iaddr,rmlexp,mptemp,lmptmp,
     $     itree,ltree,iptr,ndiv,nlevels,
     $     nboxes,iper,boxsize,rscales,centers,laddr,
     $     isrcse,itargse,iexpcse,nterms,ntj,
     $     ifpgh,pot,grad,hess,
     $     ifpghtarg,pottarg,gradtarg,hesstarg,
     $     jsort,scjsort,ifnear,timeinfo,ier,
     $     carray,ldc,
     $     mnlist1,nlist1s,list1,
     $     mnlist2,nlist2s,list2,
     $     mnlist3,nlist3s,list3,
     $     mnlist4,nlist4s,list4,
     $     nmax_plan,mploc_hexp1tmp,mploc_jexp2tmp,
     $     mploc_z0pow1,mploc_z0pow2,
     $     mpmp_z0pow1,mpmp_z0pow2,
     $     mpmp_hexp1tmp,mpmp_hexp2tmp,
     $     locloc_jexp1tmp,locloc_jexp2tmp,
     $     lmptot,box_level)
c
c   Cauchy FMM main loop using pre-computed interaction lists and
c   binomial table. Interface identical to cfmm2dmain except that
c   carray, ldc, mnlistX, nlistXs, listX are supplied by the caller
c   rather than computed internally.
c
      implicit none

      integer nd

      integer iper

      integer nsource,ntarget,nexpc
      integer ndiv,nlevels,ntj

      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      real *8 eps
      integer ier,ifnear

      real *8 sourcesort(2,nsource)

      complex *16 chargesort(nd,*)
      complex *16 dipstrsort(nd,*)

      real *8 targetsort(2,ntarget)
      complex *16 jsort(nd,0:ntj,*)

      real *8 expcsort(2,*)

      complex *16 pot(nd,*)
      complex *16 grad(nd,*)
      complex *16 hess(nd,*)

      complex *16 pottarg(nd,*)
      complex *16 gradtarg(nd,*)
      complex *16 hesstarg(nd,*)

      integer iaddr(2,nboxes),lmptmp
      real *8 rmlexp(*)
      complex *16 mptemp(lmptmp)

      real *8 timeinfo(8)
      real *8 timelev(0:200)
      real *8 centers(2,*)

      integer laddr(2,0:nlevels)
      integer nterms(0:nlevels)
      integer iptr(8),ltree
      integer itree(ltree)
      integer nboxes
      integer isrcse(2,nboxes),itargse(2,nboxes)
      integer iexpcse(2,nboxes)
      real *8 rscales(0:nlevels),boxsize(0:nlevels)

      real *8 scjsort(*)

      real *8 thresh

      integer nterms_eval(4,0:200)

c     Pre-computed list arguments
      integer ldc
      real *8 carray(0:ldc,0:ldc)
      integer mnlist1,mnlist2,mnlist3,mnlist4
      integer nlist1s(nboxes),nlist2s(nboxes)
      integer nlist3s(nboxes),nlist4s(nboxes)
      integer list1(mnlist1,nboxes),list2(mnlist2,nboxes)
      integer list3(mnlist3,nboxes),list4(mnlist4,nboxes)

c     Plan workspace arguments
      integer nmax_plan
      integer lmptot
      integer box_level(nboxes)
c     M2L step 4: per-thread scratch arrays
      complex *16 mploc_hexp1tmp(nd,0:nmax_plan,*)
      complex *16 mploc_jexp2tmp(nd,0:nmax_plan,*)
      complex *16 mploc_z0pow1(0:nmax_plan,*)
      complex *16 mploc_z0pow2(0:nmax_plan,*)
c     M2M step 3: precomputed z0pow by quadrant and per-thread scratch
      complex *16 mpmp_z0pow1(0:nmax_plan,4)
      complex *16 mpmp_z0pow2(0:nmax_plan,4)
      complex *16 mpmp_hexp1tmp(nd,0:nmax_plan,*)
      complex *16 mpmp_hexp2tmp(nd,0:nmax_plan,*)
c     L2L step 5: per-thread scratch arrays
      complex *16 locloc_jexp1tmp(nd,0:nmax_plan,*)
      complex *16 locloc_jexp2tmp(nd,0:nmax_plan,*)

c     temp variables
      integer i,j,k,l,idim,tid,iq
      integer ibox,jbox,ilev,npts
      real *8 dx,dy
      integer omp_get_thread_num
      external omp_get_thread_num
      integer nchild,nlist1,nlist2,nlist3,nlist4

      integer istart,iend,istarts,iends
      integer isstart,isend,jsstart,jsend
      integer jstart,jend
      integer istarte,iende,istartt,iendt

      integer ifprint

      integer ifhesstarg,nn
      real *8 d,time1,time2,omp_get_wtime
      real *8 tt1,tt2
      complex *16 pottmp,gradtmp,hesstmp

      double precision pi

c     ifprint is an internal information printing flag.
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, etc. if ifprint=2.

        ifprint=0

        pi = 4*atan(1.0d0)

        do i=0,nlevels
          timelev(i) = 0
        enddo

c
c     ... set the expansion coefficients to zero
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim,i,j)
      do i=1,nexpc
         do j = 0,ntj
           do idim=1,nd
             jsort(idim,j,i)=0
           enddo
         enddo
      enddo
C$OMP END PARALLEL DO
C

        do i=1,8
          timeinfo(i)=0
        enddo
c
c       ... set all multipole and local expansions to zero
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,lmptot
         rmlexp(i) = 0.0d0
      enddo
C$OMP END PARALLEL DO

c     Set scjsort
      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,i)
         do ibox = laddr(1,ilev), laddr(2,ilev)
            nchild = itree(iptr(4)+ibox-1)
            if(nchild.eq.0) then
                istart = iexpcse(1,ibox)
                iend = iexpcse(2,ibox)
                do i=istart,iend
                   scjsort(i) = rscales(ilev)
                enddo
            endif
         enddo
C$OMP END PARALLEL DO
      enddo


      if(ifprint .ge. 1)
     $   call prinf('=== STEP 1 (form mp) ====*',i,0)
        call cpu_time(time1)
C$        time1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions

      do ilev = 2,nlevels
C
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
          do ibox=laddr(1,ilev),laddr(2,ilev)
             nchild = itree(iptr(4)+ibox-1)
             istart = isrcse(1,ibox)
             iend = isrcse(2,ibox)
             npts = iend-istart+1
c              Check if current box is a leaf box
             if(nchild.eq.0.and.npts.gt.0) then
                 call l2dformmpc(nd,rscales(ilev),
     1             sourcesort(1,istart),npts,chargesort(1,istart),
     2             centers(1,ibox),nterms(ilev),
     3             rmlexp(iaddr(1,ibox)))
             endif
          enddo
C$OMP END PARALLEL DO
        endif

        if(ifdipole.eq.1.and.ifcharge.eq.0) then
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
          do ibox=laddr(1,ilev),laddr(2,ilev)
             nchild = itree(iptr(4)+ibox-1)
             istart = isrcse(1,ibox)
             iend = isrcse(2,ibox)
             npts = iend-istart+1
c              Check if current box is a leaf box
             if(nchild.eq.0.and.npts.gt.0) then
                call l2dformmpd(nd,rscales(ilev),
     1          sourcesort(1,istart),npts,dipstrsort(1,istart),
     2          centers(1,ibox),
     3          nterms(ilev),rmlexp(iaddr(1,ibox)))
             endif
          enddo
C$OMP END PARALLEL DO
        endif

        if(ifdipole.eq.1.and.ifcharge.eq.1) then
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
          do ibox=laddr(1,ilev),laddr(2,ilev)
             nchild = itree(iptr(4)+ibox-1)
             istart = isrcse(1,ibox)
             iend = isrcse(2,ibox)
             npts = iend-istart+1
c             Check if current box is a leaf box
             if(nchild.eq.0.and.npts.gt.0) then
                call l2dformmpcd(nd,rscales(ilev),
     1             sourcesort(1,istart),npts,chargesort(1,istart),
     2             dipstrsort(1,istart),
     3             centers(1,ibox),
     4             nterms(ilev),rmlexp(iaddr(1,ibox)))
             endif
          enddo
C$OMP END PARALLEL DO
        endif
      enddo


      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(1)=time2-time1

      if(ifprint.ge.1)
     $      call prinf('=== STEP 2 (form lo) ====*',i,0)
      call cpu_time(time1)
C$        time1=omp_get_wtime()
      do ilev = 2,nlevels
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,nlist4,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
           do ibox = laddr(1,ilev),laddr(2,ilev)
              npts = 0
              if(ifpghtarg.gt.0) then
                 istart = itargse(1,ibox)
                 iend = itargse(2,ibox)
                 npts = npts + iend-istart+1
              endif

              istart = iexpcse(1,ibox)
              iend = iexpcse(2,ibox)
              npts = npts + iend-istart+1

              if(ifpgh.gt.0) then
                 istart = isrcse(1,ibox)
                 iend = isrcse(2,ibox)
                 npts = npts + iend-istart+1
              endif

              if (npts .gt. 0) then
                 do i=1,nlist4s(ibox)
                    jbox = list4(i,ibox)
                    istart = isrcse(1,jbox)
                    iend = isrcse(2,jbox)
                    npts = iend-istart+1

                    call l2dformtac(nd,rscales(ilev),
     1                   sourcesort(1,istart),npts,
     2                   chargesort(1,istart),centers(1,ibox),
     3                   nterms(ilev),rmlexp(iaddr(2,ibox)))
                 enddo
              endif
           enddo
C$OMP END PARALLEL DO
        endif
        if(ifcharge.eq.0.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,nlist4,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
           do ibox = laddr(1,ilev),laddr(2,ilev)
              npts = 0
              if(ifpghtarg.gt.0) then
                 istart = itargse(1,ibox)
                 iend = itargse(2,ibox)
                 npts = npts + iend-istart+1
              endif

              istart = iexpcse(1,ibox)
              iend = iexpcse(2,ibox)
              npts = npts + iend-istart+1

              if(ifpgh.gt.0) then
                 istart = isrcse(1,ibox)
                 iend = isrcse(2,ibox)
                 npts = npts + iend-istart+1
              endif

              if (npts .gt. 0) then
                 do i=1,nlist4s(ibox)
                    jbox = list4(i,ibox)
                    istart = isrcse(1,jbox)
                    iend = isrcse(2,jbox)
                    npts = iend-istart+1

                    call l2dformtad(nd,rscales(ilev),
     1                   sourcesort(1,istart),npts,
     2                   dipstrsort(1,istart),
     3                   centers(1,ibox),nterms(ilev),
     4                   rmlexp(iaddr(2,ibox)))
                 enddo
              endif
          enddo
C$OMP END PARALLEL DO
        endif
        if(ifcharge.eq.1.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,nlist4,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
           do ibox = laddr(1,ilev),laddr(2,ilev)
              npts = 0
              if(ifpghtarg.gt.0) then
                 istart = itargse(1,ibox)
                 iend = itargse(2,ibox)
                 npts = npts + iend-istart+1
              endif

              istart = iexpcse(1,ibox)
              iend = iexpcse(2,ibox)
              npts = npts + iend-istart+1

              if(ifpgh.gt.0) then
                 istart = isrcse(1,ibox)
                 iend = isrcse(2,ibox)
                 npts = npts + iend-istart+1
              endif

              if (npts .gt. 0) then
                 do i=1,nlist4s(ibox)
                    jbox = list4(i,ibox)
                    istart = isrcse(1,jbox)
                    iend = isrcse(2,jbox)
                    npts = iend-istart+1

                    call l2dformtacd(nd,rscales(ilev),
     1                   sourcesort(1,istart),npts,
     2                   chargesort(1,istart),dipstrsort(1,istart),
     3                   centers(1,ibox),
     3                   nterms(ilev),rmlexp(iaddr(2,ibox)))
                 enddo
              endif
          enddo
C$OMP END PARALLEL DO
        endif
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1

      if(ifprint .ge. 1)
     $     call prinf('=== STEP 3 (merge mp) ====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c
      do ilev=nlevels-1,1,-1

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,i,nchild,istart,iend,npts,mptemp,tid,iq,dx,dy)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          nchild = itree(iptr(4)+ibox-1)
          tid = omp_get_thread_num()+1
          do i=1,nchild
            jbox = itree(iptr(5)+4*(ibox-1)+i-1)
            istart = isrcse(1,jbox)
            iend = isrcse(2,jbox)
            npts = iend-istart+1
            if(npts.gt.0) then
              dx = centers(1,jbox)-centers(1,ibox)
              dy = centers(2,jbox)-centers(2,ibox)
              iq = 1
              if(dx.lt.0.0d0) iq = iq+2
              if(dy.lt.0.0d0) iq = iq+1
              call l2dmpmp_work(nd,
     1           rmlexp(iaddr(1,jbox)),nterms(ilev+1),
     2           rmlexp(iaddr(1,ibox)),nterms(ilev),carray,ldc,
     3           mpmp_z0pow1(0,iq),mpmp_z0pow2(0,iq),
     4           mpmp_hexp1tmp(1,0,tid),mpmp_hexp2tmp(1,0,tid))
            endif
          enddo
        enddo
C$OMP END PARALLEL DO
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(3)=time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== Step 4 (mp to loc) ===*',i,0)
c      ... step 3, convert multipole expansions into local
c       expansions

      call cpu_time(time1)
C$    time1=omp_get_wtime()
c     Flattened over all levels 2..nlevels: one barrier instead of nlevels-1.
c     ilev is PRIVATE, looked up from box_level(ibox) precomputed at build time.
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,istart,iend,npts,i,tid,ilev)
C$OMP$SCHEDULE(DYNAMIC)
      do ibox = laddr(1,2),laddr(2,nlevels)
        ilev = box_level(ibox)
        npts = 0
        if(ifpghtarg.gt.0) then
          istart = itargse(1,ibox)
          iend = itargse(2,ibox)
          npts = npts + iend-istart+1
        endif

        istart = iexpcse(1,ibox)
        iend = iexpcse(2,ibox)
        npts = npts + iend-istart+1

        if(ifpgh.gt.0) then
          istart = isrcse(1,ibox)
          iend = isrcse(2,ibox)
          npts = npts + iend-istart+1
        endif

        if(npts.gt.0) then
          tid = omp_get_thread_num() + 1
          do i=1,nlist2s(ibox)
            jbox = list2(i,ibox)
            call l2dmploc_work(nd,rscales(ilev),
     $        centers(1,jbox),rmlexp(iaddr(1,jbox)),nterms(ilev),
     2        rscales(ilev),centers(1,ibox),rmlexp(iaddr(2,ibox)),
     3        nterms(ilev),carray,ldc,
     4        mploc_z0pow1(0,tid),mploc_z0pow2(0,tid),
     5        mploc_hexp1tmp(1,0,tid),mploc_jexp2tmp(1,0,tid))
          enddo
        endif
      enddo
C$OMP END PARALLEL DO
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(4) = time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== Step 5 (split loc) ===*',i,0)

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      do ilev = 1,nlevels-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,i,nchild,istart,iend,npts,mptemp,tid,iq,dx,dy)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          nchild = itree(iptr(4)+ibox-1)
          istart = iexpcse(1,ibox)
          iend = iexpcse(2,ibox)
          npts = iend - istart + 1


          if(ifpghtarg.gt.0) then
            istart = itargse(1,ibox)
            iend = itargse(2,ibox)
            npts = npts + iend-istart+1
          endif

          if(ifpgh.gt.0) then
            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = npts + iend-istart+1
          endif

          if(npts.gt.0) then
            tid = omp_get_thread_num()+1
            do i=1,nchild
              jbox = itree(iptr(5)+4*(ibox-1)+i-1)
              dx = centers(1,jbox)-centers(1,ibox)
              dy = centers(2,jbox)-centers(2,ibox)
              iq = 1
              if(dx.lt.0.0d0) iq = iq+2
              if(dy.lt.0.0d0) iq = iq+1
              call l2dlocloc_work(nd,
     1          rmlexp(iaddr(2,ibox)),nterms(ilev),
     2          rmlexp(iaddr(2,jbox)),nterms(ilev+1),carray,ldc,
     3          mpmp_z0pow2(0,iq),mpmp_z0pow1(0,iq),
     4          locloc_jexp1tmp(1,0,tid),locloc_jexp2tmp(1,0,tid))
            enddo
          endif
        enddo
C$OMP END PARALLEL DO
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(5) = time2-time1

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      if(ifprint.ge.1)
     $    call prinf('=== Step 6 (mp eval) ===*',i,0)

      do ilev=1,nlevels-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nlist3,istart,iend,npts,j,i,mptemp)
C$OMP$PRIVATE(jbox)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox=laddr(1,ilev),laddr(2,ilev)
          do j=iexpcse(1,ibox),iexpcse(2,ibox)
            do i=1,nlist3s(ibox)
              jbox = list3(i,ibox)
c                 shift multipole expansion directly to box
c                 for all expansion centers
              call l2dmploc(nd,rscales(ilev+1),
     $          centers(1,jbox),rmlexp(iaddr(1,jbox)),nterms(ilev+1),
     2          scjsort(j),expcsort(1,j),jsort(1,0,j),ntj,carray,ldc)
            enddo
          enddo

c              evalute multipole expansion at all targets
          istart = itargse(1,ibox)
          iend = itargse(2,ibox)
          npts = iend-istart+1

          if(ifpghtarg.eq.1) then
            do i=1,nlist3s(ibox)
              jbox = list3(i,ibox)

              call l2dmpevalp(nd,rscales(ilev+1),
     1         centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2         nterms(ilev+1),targetsort(1,istart),npts,
     3         pottarg(1,istart))
            enddo
          endif
          if(ifpghtarg.eq.2) then
            do i=1,nlist3s(ibox)
              jbox = list3(i,ibox)
              call l2dmpevalg(nd,rscales(ilev+1),
     1          centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2          nterms(ilev+1),targetsort(1,istart),npts,
     3          pottarg(1,istart),gradtarg(1,istart))
            enddo
          endif
          if(ifpghtarg.eq.3) then
            do i=1,nlist3s(ibox)
              jbox = list3(i,ibox)

              call l2dmpevalh(nd,rscales(ilev+1),
     1          centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2          nterms(ilev+1),targetsort(1,istart),npts,
     3          pottarg(1,istart),
     3          gradtarg(1,istart),hesstarg(1,istart))
            enddo
          endif


c              evalute multipole expansion at all sources
          istart = isrcse(1,ibox)
          iend = isrcse(2,ibox)
          npts = iend-istart+1


          if(ifpgh.eq.1) then
            do i=1,nlist3s(ibox)
              jbox = list3(i,ibox)
              call l2dmpevalp(nd,rscales(ilev+1),
     1           centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2           nterms(ilev+1),sourcesort(1,istart),npts,
     3           pot(1,istart))
            enddo
          endif
          if(ifpgh.eq.2) then
            do i=1,nlist3s(ibox)
              jbox = list3(i,ibox)
              call l2dmpevalg(nd,rscales(ilev+1),
     1           centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2           nterms(ilev+1),sourcesort(1,istart),npts,
     3           pot(1,istart),grad(1,istart))
            enddo
          endif
          if(ifpgh.eq.3) then
            do i=1,nlist3s(ibox)
              jbox = list3(i,ibox)
              call l2dmpevalh(nd,rscales(ilev+1),
     1           centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2           nterms(ilev+1),sourcesort(1,istart),npts,
     3           pot(1,istart),grad(1,istart),hess(1,istart))
            enddo
          endif

        enddo
C$OMP END PARALLEL DO
      enddo

 1000 continue


      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(6) = time2-time1


      if(ifprint.ge.1)
     $    call prinf('=== step 7 (eval lo) ===*',i,0)

c     ... step 7, evaluate all local expansions
      call cpu_time(time1)
C$    time1=omp_get_wtime()
      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,mptemp,istart,iend,i,npts)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          nchild = itree(iptr(4)+ibox-1)
          if(nchild.eq.0) then
            istart = iexpcse(1,ibox)
            iend = iexpcse(2,ibox)
            do i=istart,iend
              call l2dlocloc(nd,rscales(ilev),
     $          centers(1,ibox),
     1          rmlexp(iaddr(2,ibox)),nterms(ilev),scjsort(i),
     2          expcsort(1,i),jsort(1,0,i),ntj,carray,ldc)
            enddo
c
cc               evaluate local expansion
c                at targets
            istart = itargse(1,ibox)
            iend = itargse(2,ibox)
            npts = iend-istart + 1
            if(ifpghtarg.eq.1) then
              call l2dtaevalp(nd,rscales(ilev),
     1              centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2              nterms(ilev),targetsort(1,istart),npts,
     3              pottarg(1,istart))
            endif
            if(ifpghtarg.eq.2) then
              call l2dtaevalg(nd,rscales(ilev),
     1          centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2          nterms(ilev),targetsort(1,istart),npts,
     3          pottarg(1,istart),gradtarg(1,istart))
            endif
            if(ifpghtarg.eq.3) then
              call l2dtaevalh(nd,rscales(ilev),
     1          centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2          nterms(ilev),targetsort(1,istart),npts,
     3          pottarg(1,istart),gradtarg(1,istart),
     4          hesstarg(1,istart))
            endif

c
cc                evaluate local expansion at sources

            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = iend-istart+1
            if(ifpgh.eq.1) then
              call l2dtaevalp(nd,rscales(ilev),
     1           centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2           nterms(ilev),sourcesort(1,istart),npts,
     3           pot(1,istart))
            endif
            if(ifpgh.eq.2) then
              call l2dtaevalg(nd,rscales(ilev),
     1           centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2           nterms(ilev),sourcesort(1,istart),npts,
     3           pot(1,istart),grad(1,istart))
            endif
            if(ifpgh.eq.3) then
              call l2dtaevalh(nd,rscales(ilev),
     1           centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2           nterms(ilev),sourcesort(1,istart),npts,
     3           pot(1,istart),grad(1,istart),hess(1,istart))
            endif
          endif
        enddo
C$OMP END PARALLEL DO
      enddo

      call cpu_time(time2)
C$    time2 = omp_get_wtime()
      timeinfo(7) = time2 - time1

      if(ifprint .ge. 1)
     $     call prinf('=== STEP 8 (direct) =====*',i,0)

c
cc     set threshold for ignoring interactions with
c      |r| < thresh
c
      thresh = boxsize(0)*2.0d0**(-51)

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      if(ifnear.eq.0) goto 1233

      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,istartt,iendt,i,jstart,jend,istarte,iende)
C$OMP$PRIVATE(nlist1,istarts,iends)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = laddr(1,ilev),laddr(2,ilev)

            istartt = itargse(1,ibox)
            iendt = itargse(2,ibox)


            istarte = iexpcse(1,ibox)
            iende = iexpcse(2,ibox)

            istarts = isrcse(1,ibox)
            iends = isrcse(2,ibox)

            do i =1,nlist1s(ibox)
               jbox = list1(i,ibox)

               jstart = isrcse(1,jbox)
               jend = isrcse(2,jbox)

               call cfmm2dexpc_direct(nd,jstart,jend,istarte,
     1         iende,rscales,nlevels,
     2         sourcesort,ifcharge,chargesort,ifdipole,dipstrsort,
     3         expcsort,jsort,scjsort,ntj)


               call cfmm2dpart_direct(nd,jstart,jend,istartt,
     1         iendt,sourcesort,ifcharge,chargesort,ifdipole,
     2         dipstrsort,targetsort,ifpghtarg,pottarg,
     3         gradtarg,hesstarg,thresh)

               call cfmm2dpart_direct(nd,jstart,jend,istarts,iends,
     1         sourcesort,ifcharge,chargesort,ifdipole,
     2         dipstrsort,sourcesort,ifpgh,pot,grad,hess,
     3         thresh)
            enddo
         enddo
C$OMP END PARALLEL DO
      enddo
 1233 continue
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(8) = time2-time1
      if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,8)
      d = 0
      do i = 1,8
         d = d + timeinfo(i)
      enddo

      if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)
      if(ifprint.ge.1) call prin2('timlev=*',timelev,nlevels+1)

      return
      end
