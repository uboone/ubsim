# ======================================================================
#  larsoft main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larsoft/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================

cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

project(ubsim VERSION 09.65.03 LANGUAGES CXX)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

list(APPEND CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules $ENV{ART_DIR}/Modules)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product( larcorealg )
find_ups_product( larcoreobj )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( ubobj )
find_ups_boost( )
find_ups_root( )
find_ups_product( canvas )
find_ups_product( ubevt )
find_ups_product( ubcore )
find_ups_product( larsim )
find_ups_product( larevt )
find_ups_product( lardata )
find_ups_product( art )
find_ups_product(art_root_io)
find_ups_product( larcore )
find_ups_product( nutools )
find_ups_product( nug4 )
find_ups_product( nugen )
find_ups_product( nurandom )
find_ups_product( ifdh_art )
find_ups_product( genie )

# WireCell required
find_ups_product( jsoncpp )
find_ups_product( gojsonnet )
find_ups_product( spdlog )
find_ups_product( wirecell )
find_ups_product( eigen )
find_ups_product( larwirecell )

cet_find_library( WIRECELL_SIGPROC_LIB NAMES WireCellSigProc PATHS ENV WIRECELL_LIB NO_DEFAULT_PATH )
cet_find_library( WIRECELL_GEN_LIB NAMES WireCellGen PATHS ENV WIRECELL_LIB NO_DEFAULT_PATH )
cet_find_library( WIRECELL_IFACE_LIB NAMES WireCellIface PATHS ENV WIRECELL_LIB NO_DEFAULT_PATH )
cet_find_library( WIRECELL_UTIL_LIB NAMES WireCellUtil PATHS ENV WIRECELL_LIB NO_DEFAULT_PATH )
cet_find_library( WIRECELL_APPS_LIB NAMES WireCellApps PATHS ENV WIRECELL_LIB NO_DEFAULT_PATH )
cet_find_library( WIRECELL_AUX_LIB NAMES WireCellAux PATHS ENV WIRECELL_LIB NO_DEFAULT_PATH )
set(WIRECELL_LIBS ${WIRECELL_APPS_LIB} ${WIRECELL_SIGPROC_LIB} ${WIRECELL_IFACE_LIB} ${WIRECELL_UTIL_LIB} ${WIRECELL_GEN_LIB} ${WIRECELL_AUX_LIB})

cet_find_library( JSONCPP NAMES jsoncpp PATHS ENV JSONCPP_LIB NO_DEFAULT_PATH )

# geant4reweight required
find_ups_product( geant4reweight )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(ubsim)

# Unit tests.
add_subdirectory(test)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
