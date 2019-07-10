install_fhicl()
FILE(GLOB fcl_files *.fcl)
install_source( EXTRAS ${fcl_files} )

add_subdirectory(flux)
add_subdirectory(genie)
add_subdirectory(xs)
add_subdirectory(reint)
add_subdirectory(lee_weights)
add_subdirectory(splines)

