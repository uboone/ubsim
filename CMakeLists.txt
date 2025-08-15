# ======================================================================
#  ubsim main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/lardataobj/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules REQUIRED)
project(ubsim VERSION 10.10.00 LANGUAGES CXX)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
              $<$<CXX_COMPILER_ID:Clang>:-Wno-c99-extensions>
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

find_package( ubevt REQUIRED EXPORT )
find_package( ubcore REQUIRED EXPORT )
find_package( larsim REQUIRED EXPORT )
find_package(art_root_io REQUIRED EXPORT )
find_package( nurandom REQUIRED EXPORT )
find_package(nugen REQUIRED)
find_package(nutools REQUIRED)
find_package(ppfx REQUIRED)
find_package( ifdh_art REQUIRED EXPORT )
find_package( larwirecell REQUIRED EXPORT )
find_package( larcore  REQUIRED EXPORT )
find_package( larcorealg  REQUIRED EXPORT )
find_package( geant4reweight REQUIRED EXPORT )
find_package(LibXml2 REQUIRED EXPORT)
find_package(log4cpp REQUIRED EXPORT)

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(ubsim)

# Unit tests.
add_subdirectory(test)

# ups - table and config files

# packaging utility
cet_cmake_config()
