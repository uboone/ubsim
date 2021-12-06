
art_make(
	  LIB_LIBRARIES
			#messagefacility::MF_MessageLogger
			#
			#cetlib::cetlib cetlib_except

		        art::Framework_Core
		        fhiclcpp::fhiclcpp
        )


install_headers()
install_fhicl()
install_source()

