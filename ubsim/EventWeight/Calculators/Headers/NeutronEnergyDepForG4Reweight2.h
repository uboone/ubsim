
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

std::vector< std::pair<double, int> > NeutronEnergyDepForG4Reweight2(G4ReweightTraj * theTraj, double res, double mass, bool isElastic){

  std::vector< std::pair<double, int> > result;

  //First slice position
//  double sliceEdge = res;
//  double lastPos = 0.;
//  double nextPos = 0.;
//  double px,py,pz;
  int interactInSlice = 0;

  //Get total distance traveled in z
//  double totalDeltaZ = 0.;
  double disp = 0.;
//  double oldDisp = 0.;
//  int crossedSlices = 0;

  int currentSlice = 0;
  int oldSlice = 0;

double stepEnergyLoss = 0;

  double sliceEnergy = theTraj->GetEnergy();
  size_t nSteps = theTraj->GetNSteps();
  for(size_t is = 0; is < nSteps; ++is){
    auto theStep = theTraj->GetStep(is);

    // Energy at start of step
    double stepBeginEnergy = sqrt( theTraj->GetStep(is)->GetPreStepPx()*theTraj->GetStep(is)->GetPreStepPx()+theTraj->GetStep(is)->GetPreStepPy()*theTraj->GetStep(is)->GetPreStepPy()+theTraj->GetStep(is)->GetPreStepPz()*theTraj->GetStep(is)->GetPreStepPz()+mass*mass);

    // Energy at end of step
    double stepEndEnergy = sqrt( theTraj->GetStep(is)->GetPostStepPx()*theTraj->GetStep(is)->GetPostStepPx()+theTraj->GetStep(is)->GetPostStepPy()*theTraj->GetStep(is)->GetPostStepPy()+theTraj->GetStep(is)->GetPostStepPz()*theTraj->GetStep(is)->GetPostStepPz()+mass*mass);

         stepEnergyLoss = stepBeginEnergy - stepEndEnergy;

    disp += theStep->GetStepLength();
    currentSlice = floor(disp/res);

    std::string theProc = theStep->GetStepChosenProc();

    //Check to see if in a new slice and it's not the end
    if( oldSlice != currentSlice && is < nSteps - 1){

      //Save Interaction info of the prev slice
      //and reset
      result.push_back( std::make_pair(sliceEnergy, interactInSlice) );
      interactInSlice = 0;

      //Update the energy
      //sliceEnergy = sliceEnergy - res*BetheBloch(sliceEnergy, mass);
      //sliceEnergy = sliceEnergy - theStep->GetStepLength()/res*stepEnergyLoss;
      sliceEnergy = sliceEnergy - res/theStep->GetStepLength()*stepEnergyLoss;

      if( sliceEnergy - mass < 0.){
        //std::cout << "Warning! Negative energy " << sliceEnergy - mass << std::endl;
        //std::cout << "Crossed " << oldSlice - currentSlice << std::endl;
        sliceEnergy = 0.0001;
      }
      //If it's more than 1 slice, add in non-interacting slices
      for(int ic = 1; ic < abs( oldSlice - currentSlice ); ++ic){
        //std::cout << ic << std::endl;

        result.push_back( std::make_pair(sliceEnergy, 0) );

        //Update the energy again
        //sliceEnergy = sliceEnergy - res*BetheBloch(sliceEnergy, mass);
        //sliceEnergy = sliceEnergy - theStep->GetStepLength()/res*stepEnergyLoss;
        sliceEnergy = sliceEnergy - res/theStep->GetStepLength()*stepEnergyLoss;
        if( sliceEnergy - mass < 0.){
          //std::cout << "Warning! Negative energy " << sliceEnergy - mass << std::endl;
          //std::cout << "Crossed " << oldSlice - currentSlice << std::endl;
          sliceEnergy = 0.0001;
        }
      }

      if ((!isElastic && theProc.find(std::string("Inelastic")) != std::string::npos) || (isElastic && theProc.find(std::string("hadElastic")) != std::string::npos)) {
          //std::cout << "found! " << theProc << '\n';
          interactInSlice = 1;
      }
    }
    //It's crossed a slice and it's the last step. Save both info
    else if( oldSlice != currentSlice && is == nSteps - 1 ){
      result.push_back( std::make_pair(sliceEnergy, interactInSlice) );
      interactInSlice = 0;

      //Update the energy
        ///sliceEnergy = sliceEnergy - theStep->GetStepLength()/res*stepEnergyLoss;
        sliceEnergy = sliceEnergy - res/theStep->GetStepLength()*stepEnergyLoss;
      //sliceEnergy = sliceEnergy - res*BetheBloch(sliceEnergy, mass);
      if( sliceEnergy - mass < 0.){
        //std::cout << "Warning! Negative energy " << sliceEnergy - mass << std::endl;
        //std::cout << "Crossed " << oldSlice - currentSlice << std::endl;
        sliceEnergy = 0.0001;
      }
      //If it's more than 1 slice, add in non-interacting slices
      for(int ic = 1; ic < abs( oldSlice - currentSlice ); ++ic){
        //std::cout << ic << std::endl;

        result.push_back( std::make_pair(sliceEnergy, 0) );

        //Update the energy again
        //sliceEnergy = sliceEnergy - theStep->GetStepLength()/res*stepEnergyLoss;
        sliceEnergy = sliceEnergy - res/theStep->GetStepLength()*stepEnergyLoss;
        //sliceEnergy = sliceEnergy - res*BetheBloch(sliceEnergy, mass);
        if( sliceEnergy - mass < 0.){
          //std::cout << "Warning! Negative energy " << sliceEnergy - mass << std::endl;
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
    }
    //It's the end, so just save this last info
    else if( oldSlice == currentSlice && is == nSteps - 1 ){
      if ((!isElastic && theProc.find(std::string("Inelastic")) != std::string::npos) || (isElastic && theProc.find(std::string("hadElastic")) != std::string::npos)) {
              // std::cout << "found! " << theProc << '\n';
              interactInSlice = 1;
          }
      result.push_back( std::make_pair(sliceEnergy, interactInSlice) );
    }
    //Same slice, not the end. Check for interactions
    else{
      if ((!isElastic && theProc.find(std::string("Inelastic")) != std::string::npos) || (isElastic && theProc.find(std::string("hadElastic")) != std::string::npos)) {
          //std::cout << "found! " << theProc << '\n';
          interactInSlice = 1;
      }
    }

    //Update oldslice
    oldSlice = currentSlice;
  }

  return result;
}
