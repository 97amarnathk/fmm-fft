      subroutine myffti(len, wkkeep, wktemp, szkeep, sztemp)
C
C     Initialize work arrays for a sequential one-dimensional FFT routine.
C
C     integer len:  length of transform
C
C     double precision wkkeep:  work array to be initialized,
C                               not to be modified between calls to myfft
C
C     double precision wktemp:  temporary work array,
C                               may be modified between calls to myfft
C
C     integer szkeep:  size of wkkeep (at least 1)
C                      This depends on len and the specific FFT routine;
C                      for dcffti/dcfftf, it is 4*len+15.
C
C     integer sztemp:  size of wktemp (at least 1)
C                      This depends on len and the specific FFT routine;
C                      for dcffti/dcfftf, it is 1.
C
      implicit none
      integer len, szkeep, sztemp
      double precision wkkeep(szkeep), wktemp(sztemp)

      call dcffti(len, wkkeep)

      return
      end


      subroutine myfft(len, x, wkkeep, wktemp, szkeep, sztemp)
C
C     Perform forward one-dimensional FFT using sequential library routine.
C
C     integer len:  length of transform
C
C     double complex x:  input/output vector of length len
C
C     double precision wkkeep:  work array initialized in myffti,
C                               not to be modified between calls to myfft
C
C     double precision wktemp:  temporary work array,
C                               may be modified between calls to myfft
C
C     integer szkeep:  size of wkkeep (at least 1)
C                      This depends on len and the specific FFT routine;
C                      for dcffti/dcfftf, it is 4*len+15.
C
C     integer sztemp:  size of wktemp (at least 1)
C                      This depends on len and the specific FFT routine;
C                      for dcffti/dcfftf, it is 1.
C
      implicit none
      integer len, szkeep, sztemp
      double precision wkkeep(szkeep), wktemp(sztemp)
      double complex x(len)

      call dcfftf(len, x, wkkeep)

      return
      end
