# for classes_def.xml!!

cet_make( LIBRARIES ${ROOT_BASIC_LIB_LIST}
	  	    ${ROOT_EG}
		    NO_DICTIONARY
)

art_dictionary( DICTIONARY_LIBRARIES ubooneobj_Optical )

install_headers()
install_fhicl()
install_source()
