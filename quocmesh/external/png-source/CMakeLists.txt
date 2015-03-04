# Silence all GCC warnings triggered by libpng.
IF ( CMAKE_COMPILER_IS_GNUCXX )
  ADD_DEFINITIONS ( -w )
  REMOVE_DEFINITIONS ( -pedantic )
ENDIF ( CMAKE_COMPILER_IS_GNUCXX )

SET ( HDRS
  png.h
  pngconf.h
)

SET ( SRCS
  example.c
  png.c
  pngerror.c
  pnggccrd.c
  pngget.c
  pngmem.c
  pngpread.c
  pngread.c
  pngrio.c
  pngrtran.c
  pngrutil.c
  pngset.c
  pngtrans.c
  pngvcrd.c
  pngwio.c
  pngwrite.c
  pngwtran.c
  pngwutil.c
)

SET ( EXTERNAL_PNG_INCLUDE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}
  CACHE INTERNAL ""
)

ADD_LIBRARY ( png STATIC ${SRCS} ${HDRS} )

SET ( EXTERNAL_PNG_LIBRARIES png CACHE INTERNAL "" )