simple_plugin(SparseRawDigitTest "module"
              larcore_Geometry
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_BASIC}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
              ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
              art_Persistency_Common canvas_Persistency_Common
              art_Persistency_Provenance canvas_Persistency_Provenance
              art_Utilities canvas_Utilities
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${FHICLCPP}
              ${CETLIB}
              ${PQ}
              BASENAME_ONLY
             )

install_headers()
install_fhicl()
install_source()
