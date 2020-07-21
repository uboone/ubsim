
//don't need this right now
/*
 double BetheBloch(double energy, double mass){

  //Need to make this configurable? Or delete...
  double K = .307075;
  double rho = 1.390;
  double Z = 18;
  double A = 40;
  double I = 188E-6;
  double me = .511;
  //Need to make sure this is total energy, not KE
  double gamma = energy/mass;
  double beta = sqrt( 1. - (1. / (gamma*gamma)) );
  double Tmax = 2 * me * beta*beta * gamma*gamma;

  double first = K * (Z/A) * rho / (beta*beta);
  double second = .5 * log(Tmax*Tmax/(I*I)) - beta*beta;

  double dEdX = first*second;
  return dEdX;
}
*/



std::vector< std::pair<double, int> > NeutronEnergyDep(G4ReweightTraj * theTraj, double res, double mass, bool isElastic ){

  std::vector< std::pair<double, int> > result;


 int interactInSlice = 0;

  double disp = 0.;

  int currentSlice = 0;
  int oldSlice = 0;

//get initial energy
double sliceEnergy = theTraj->GetEnergy();
double stepEnergyLoss=0;

double stepBeginEnergy;
double stepEndEnergy;

int slicesInStep=0;
int calcSlicesInStep=0;
double sliceEnergyLoss=0; 


 size_t nSteps = theTraj->GetNSteps();
  for(size_t is = 0; is < nSteps; ++is){
    auto theStep = theTraj->GetStep(is);
  
  disp += theStep->GetStepLength();


//	steplength = theStep->GetStepLength();


//calculate energy loss over entrie step
if(is < nSteps-1){	    

	stepBeginEnergy = sqrt( theTraj->GetStep(is)->GetPreStepPx()*theTraj->GetStep(is)->GetPreStepPx()+theTraj->GetStep(is)->GetPreStepPy()*theTraj->GetStep(is)->GetPreStepPy()+theTraj->GetStep(is)->GetPreStepPz()*theTraj->GetStep(is)->GetPreStepPz()+mass*mass);

	stepEndEnergy = sqrt( theTraj->GetStep(is+1)->GetPreStepPx()*theTraj->GetStep(is+1)->GetPreStepPx()+theTraj->GetStep(is+1)->GetPreStepPy()*theTraj->GetStep(is+1)->GetPreStepPy()+theTraj->GetStep(is+1)->GetPreStepPz()*theTraj->GetStep(is+1)->GetPreStepPz()+mass*mass);


	stepEnergyLoss = stepBeginEnergy - stepEndEnergy;
//std::cout << "Slice energy at begin step: "  << sliceEnergy << std::endl;
//std::cout << "Start of step:  " <<  stepBeginEnergy << std::endl;
//std::cout << "End of step:  " << stepEndEnergy << std::endl;
//std::cout << "Energy loss:  " << stepEnergyLoss << std::endl;


}
else
	stepEnergyLoss=0;


//Count num of slices in the step

currentSlice = floor(disp/res);


    std::string theProc = theStep->GetStepChosenProc();

calcSlicesInStep= abs( oldSlice - currentSlice );

//assume energy loss in each slice is equal to total energy loss / num of slices
if(calcSlicesInStep >0) sliceEnergyLoss =  stepEnergyLoss / calcSlicesInStep;
else sliceEnergyLoss=0;
  
  //Check to see if in a new slice and it's not the end
    if( oldSlice != currentSlice && is < nSteps - 1){


      //Save Interaction info of the prev slice
      //and reset
     result.push_back( std::make_pair(sliceEnergy, interactInSlice) );

      interactInSlice = 0;

//find new way to compute energy loss for slice


    
sliceEnergy = sliceEnergy - sliceEnergyLoss;
	slicesInStep++;	      
if( sliceEnergy - mass < 0.){
       std::cout << "Warning! Negative energy1  " << sliceEnergy - mass << std::endl;
        //std::cout << "Crossed " << oldSlice - currentSlice << std::endl;
    //    sliceEnergy = 0.0001;
      }
      


//calcSlicesInStep= abs( oldSlice - currentSlice );


//If it's more than 1 slice, add in non-interacting slices
      for(int ic = 1; ic < abs( oldSlice - currentSlice ); ++ic){

        result.push_back( std::make_pair(sliceEnergy, 0) );

	slicesInStep++;
        //Update the energy again
        sliceEnergy = sliceEnergy - sliceEnergyLoss;
 
       if( sliceEnergy - mass < 0.){
       std::cout << "Warning! Negative energy2  " << sliceEnergy - mass << std::endl;
          //std::cout << "Crossed " << oldSlice - currentSlice << std::endl;
          sliceEnergy = 0.0001;
        }
      }//for(int ic=1;ic<abs(oldSlice-currentSlice); ++ic)

      if ((!isElastic && theProc.find(std::string("Inelastic")) != std::string::npos) || (isElastic && theProc.find(std::string("hadElastic")) != std::string::npos)) {
          // std::cout << "found! " << theProc << '\n';
     //     interactInSlice = 1;
      }
    }// if(oldSlice != currentSile && is < nSteps-1)





    //It's crossed a slice and it's the last step. Save both info
    else if( oldSlice != currentSlice && is == nSteps - 1 ){
      result.push_back( std::make_pair(sliceEnergy, interactInSlice) );
      interactInSlice = 0;
//		std::cout << "Slice Energy: " << sliceEnergy << std::endl;

	slicesInStep++;
      //Update the energy
      sliceEnergy = sliceEnergy - sliceEnergyLoss;
      if( sliceEnergy - mass < 0.){
        std::cout << "Warning! Negative energy3  " << sliceEnergy - mass << std::endl;
        //std::cout << "Crossed " << oldSlice - currentSlice << std::endl;
        sliceEnergy = 0.0001;
      }
      //If it's more than 1 slice, add in non-interacting slices
      for(int ic = 1; ic < abs( oldSlice - currentSlice ); ++ic){
        //std::cout << ic << std::endl;

        result.push_back( std::make_pair(sliceEnergy, 0) );
//			std::cout << "Slice Energy: " << sliceEnergy << std::endl;

	slicesInStep++;
        //Update the energy again
        sliceEnergy = sliceEnergy - sliceEnergyLoss;
        if( sliceEnergy - mass < 0.){
          std::cout << "Warning! Negative energy4  " << sliceEnergy - mass << std::endl;
          //std::cout << "Crossed " << oldSlice - currentSlice << std::endl;
          sliceEnergy = 0.0001;
        }
      }

      //Save the last slice
      if ((!isElastic && theProc.find(std::string("Inelastic")) != std::string::npos) || (isElastic && theProc.find(std::string("hadElastic")) != std::string::npos)) {
          // std::cout << "found! " << theProc << '\n';
          interactInSlice = 1;
      }
      result.push_back( std::make_pair(sliceEnergy, interactInSlice) );
//			std::cout << "Slice Energy: " << sliceEnergy << std::endl;
  
  }//if(oldSlice != currentSlice && is == nSteps-1)





    //It's the end, so just save this last info
    else if( oldSlice == currentSlice && is == nSteps - 1 ){
      if ((!isElastic && theProc.find(std::string("Inelastic")) != std::string::npos) || (isElastic && theProc.find(std::string("hadElastic")) != std::string::npos)) {
              // std::cout << "found! " << theProc << '\n';
              interactInSlice = 1;
          }
     result.push_back( std::make_pair(sliceEnergy, interactInSlice) );
//		std::cout << "Slice Energy: " << sliceEnergy << std::endl;

    }//else if(oldSlice == currentSlice && is==nSteps-1)



    //Same slice, not the end. Check for interactions
    else{
      if ((!isElastic && theProc.find(std::string("Inelastic")) != std::string::npos) || (isElastic && theProc.find(std::string("hadElastic")) != std::string::npos)) {
          // std::cout << "found! " << theProc << '\n';
          interactInSlice = 1;
      }
    }//else

    //Update oldslice
    oldSlice = currentSlice;

if(calcSlicesInStep != slicesInStep){

std::cout << "Est slices in this step: " << calcSlicesInStep << std::endl;
std::cout <<"Slices in this step: " <<  slicesInStep << std::endl;


}

slicesInStep=0;




  }//loop over nSteps

  return result;
}
