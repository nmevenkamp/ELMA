SET ( HDRS
  _kiss_fft_guts.h
  kiss_fft.h
  kissfft.hh
  kfc.h
  kiss_fftnd.h
  kiss_fftndr.h
  kiss_fftr.h
)

SET ( SRCS
  kiss_fft.c
  kfc.c
  kiss_fastfir.c
  kiss_fftnd.c
  kiss_fftndr.c
  kiss_fftr.c
)

SET ( EXTERNAL_KISSFFT_INCLUDE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}
  CACHE INTERNAL ""
)

ADD_LIBRARY ( kissfft STATIC ${SRCS} ${HDRS} )

SET ( EXTERNAL_KISSFFT_LIBRARIES kissfft CACHE INTERNAL "" )