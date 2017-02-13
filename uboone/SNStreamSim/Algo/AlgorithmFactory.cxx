#ifndef COMPRESS_ALGORITHMFACTORY_CXX
#define COMPRESS_ALGORITHMFACTORY_CXX

#include "AlgorithmFactory.h"

namespace compress {

  AlgoDefault::AlgoDefault() {

    std::cout << "setting up algorithm..." << std::endl;
    
    this->SetCompressThresh(30,30,30);
    this->SetPolarity(1,1,1);
    this->SetMaxADC(4095);
    this->SetUVYplaneBuffer(7,7,7,7,7,7);
    this->SetBlockSize(64);
    this->SetBaselineThresh(2.0);
    this->SetVarianceThresh(2.0);
    this->SetDebug(false);

  }

}

#endif
