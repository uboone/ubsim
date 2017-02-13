#ifndef COMPRESS_ALGORITHMFACTORY_H
#define COMPRESS_ALGORITHMFACTORY_H

#include <string>
#include <cstdlib>
#include "uboone/SNStreamSim/Fmwk/CompressionAlgoBase.h"
#include "CompressionAlgosncompress.h"

namespace compress {


  class AlgoDefault : public CompressionAlgosncompress {

  public:
    AlgoDefault();
    
    ~AlgoDefault(){}

  };

}

#endif
