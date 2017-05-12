#include "uboone/CRT/CRTAuxFunctions.hh"

int crt::auxfunctions::getFEBN(uint64_t febid){

  int FEBN = (febid  - 413240800282 + 26)%100; //FEB ID                                                                                                                 
  return FEBN;
}


void crt::auxfunctions::FillPos(std::string filePos,   std::map <int, std::vector<double> >& sensor_pos){ //key = FEB*100+ch                                                   

  std::cout<<"Reading SiPMs position information"<<std::endl;

  std::ifstream in;
  in.open(filePos.c_str());

  std::cout <<"File open? "<<in.is_open()<<std::endl;

  if(in.is_open()){
    std::cout<<"File open: "<<filePos.c_str()<<std::endl;
  }

  // getchar();

  while (!in.eof()){

    double x,y,z;
    int id, pl, subpl, ori;
    std::vector<double> pos;

    in>>id>>x>>y>>z>>pl>>subpl>>ori;

    pos.push_back(x);
    pos.push_back(y);
    pos.push_back(z);
    pos.push_back(pl);
    pos.push_back(subpl);
    pos.push_back(ori);

    sensor_pos[id] = pos;


  }

}

std::vector<double> crt::auxfunctions::getPos(int ID,  std::map <int, std::vector<double> >& sensor_pos){

  auto it2 = sensor_pos.find(ID);
  std::vector<double> pos = (*it2).second;

  return pos;
  
}



void crt::auxfunctions::FillFEBDel(std::string fileFEBDel,   std::map <int,double> & FEBDel){ //key = FEB                                                                                       

  std::cout<<"Reading FEB delay information"<<std::endl;
 
  std::ifstream in;
  in.open(fileFEBDel.c_str());

  std::cout <<"File open? "<<in.is_open()<<std::endl;

  if(in.is_open()){
    std::cout<<"File open: "<<fileFEBDel.c_str()<<std::endl;
  }

  // getchar();                                                                                                                                                                                                     

  while (!in.eof()){

    double del;
    int id;

    in>>id>>del;
    FEBDel[id] = del;

  }

}

double crt::auxfunctions::getFEBDel(int ID,  std::map <int,double >& FEBDel){

  auto it2 = FEBDel.find(ID);
  double del = (*it2).second;

  return del;

}


//
void crt::auxfunctions::FillGain(std::string fileGain,   std::map <int, std::pair<double,double> >& sensor_gain){ //key = FEB*100+ch            

  std::cout<<"Reading SiPMs gain information"<<std::endl;
  

  std::ifstream in;
  in.open(fileGain.c_str());

  std::cout <<"File open? "<<in.is_open()<<std::endl;

  if(in.is_open()){
    std::cout<<"File open: "<<fileGain.c_str()<<std::endl;
  }

  // getchar();                                                                                                                                                    

  while (!in.eof()){

    double g, g_err;
    int id;
    std::pair<double,double> gain;

    in>>id>>g>>g_err;
    gain = std::make_pair(g,g_err);
    sensor_gain[id] = gain;

                                                                                                                                                                     
  }

}


std::pair<double,double> crt::auxfunctions::getGain(int ID,  std::map <int, std::pair<double,double> >& sensor_gain){

  auto it2 = sensor_gain.find(ID);
  std::pair<double,double> gain = (*it2).second;
  return gain;
}

//
std::vector<double> crt::auxfunctions::inter_X(double S1, std::vector<double>& pS1, double S2, std::vector<double>& pS2){

  std::vector<double> interpos;

  //  std::cout<<"  pos1 x:  "<<pS1[0]<<" y: "<<pS1[1]<<" z: "<<pS1[2]<<std::endl;
  //std::cout<<"  pos2 x:  "<<pS2[0]<<" y: "<<pS2[1]<<" z: "<<pS2[2]<<std::endl;
  //std::cout<<" "<<std::endl;

  if(S1<0) S1 = 0.000000001;
  if(S2<0) S2 = 0.000000001;
  
  //double dy = 0.344677; //pol1
  //double d0y = -1.92045;
  
  /*  double p0y = -3.14878; //pol3
  double p1y =  1.15249;
  double p2y =  -0.154124;
  double p3y = 0.00892841;
  */
  /*
  double x =  ((S1*0 + S2*L)/(S1+S2)); //barycenter
  x =  x + dy * x + d0y; //barycenter Norm pol1
   
  if(x<0) x=0;
  if(x>L) x=L;
*/  

  double x =  ((S1*pS1[0] + S2*pS2[0])/(S1+S2)); //barycenter
  // x =  x + dy * x + d0y; //barycenter Norm pol
  //std::cout<<"x*: "<<x<<std::endl;
  //if(x<pS1[0]) x=pS1[0];
  //if(x>pS2[0]) x=pS2[0];
  //std::cout<<"x: "<<x<<std::endl;

  double y =  ((S1*pS1[1] + S2*pS2[1])/(S1+S2)); //barycenter
  //y =  y + dy * y + d0y; //barycenter Norm pol1
  //std::cout<<"y*: "<<y<<std::endl;
  //if(y<pS1[0]) y=pS1[1];
  //if(y>pS2[0]) y=pS2[1];
  //std::cout<<"y: "<<y<<std::endl;

  double z =  ((S1*pS1[2] + S2*pS2[2])/(S1+S2)); //barycenter
  //z =  z + dy * z + d0y; //barycenter Norm pol1
  // std::cout<<"z*: "<<z<<std::endl;
  //if(z<pS1[2]) z=pS1[2];
  //if(z>pS2[2]) z=pS2[2];
  //std::cout<<"z: "<<z<<std::endl;

  double pl = (pS1[3] + pS2[3])/2;
  double ly = (pS1[4] + pS2[4])/2;
  double ori = (pS1[5] + pS2[5])/2;
  double inter = 0;

  if(pS1[0] != pS2[0]) inter=1; //interpol x
  if(pS1[1] != pS2[1]) inter=2; //interpol y
  if(pS1[2] != pS2[2]) inter=3; //interpol z

  interpos.push_back(x);
  interpos.push_back(y);
  interpos.push_back(z);
  interpos.push_back(pl);
  interpos.push_back(ly);
  interpos.push_back(ori);
  interpos.push_back(inter);

  return interpos;
}

double crt::auxfunctions::inter_X_error(double S1, double S2, double L){

  if(S1<0) S1 = 0.000000001;
  if(S2<0) S2 = 0.000000001;
  
  double dy = 0.344677; //pol1
  double d0y = -1.92045;
  /*  
  double p0y = -3.14878; //pol3
  double p1y =  1.15249;
  double p2y =  -0.154124;
  double p3y = 0.00892841;
  */
  double x =  ((S1*0 + S2*L)/(S1+S2)); //barycenter
  x =  x + dy * x + d0y; //barycenter Norm pol1

  double x_error =  1.92380e+00+1.47186e-02*x-5.29446e-03*x*x;

  //std::cout<<"pos: "<<x<<std::endl;
  //std::cout<<"error: "<<x_error<<std::endl;
  
  if(x_error<1.45) x_error=1.45;
  if(x_error>1.99) x_error=1.99;

  return x_error;
}


double crt::auxfunctions::getTcorr(std::vector<double>& inpos1, std::vector<double>& inpos2 ,double T){                                                                         
  
  double L = 0;                                                                                                                                              
  if(inpos1[3]==0 || inpos1[3]==3){//bottom&top (XZ)                                                                                                         
    
    if(inpos1[6]==1){//correction in z                                                                                                                       
      L = abs(inpos1[2]-inpos2[2]);                                                                                                                          
    }                                                                                                                                                        
    
    if(inpos1[6]==3){//correction in x                                                                                                                       
      L = abs(inpos1[0]-inpos2[0]);                                                                                                                          
    }                                                                                                                                                        
  }                                                                                                                                                          
  
  if(inpos1[3]==1 || inpos1[3]==2){//FT&pipe (YZ)                                                                                                           
    
    if(inpos1[6]==2){//correction in z                                                                                                                       
      L = abs(inpos1[2]-inpos2[2]);                                                                                                                          
    }                                                                                                                                                        
    
    if(inpos1[6]==3){//correction in y                                                                                                                       
      L = abs(inpos1[1]-inpos2[1]);                                                                                                                          
    }                                                                                                                                                        
    
  }                                                                                                                                                          
  
  double Tcorr = T - (L*6.2/100);                                                                                                                            
  return Tcorr;                                                                                                                                              
}
