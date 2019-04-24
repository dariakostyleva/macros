{
00191   if (GetVerboseLevel()>1) G4cout << "G4GeneralPhaseSpaceDecay::TwoBodyDecayIt()"<<G4endl;
00192   
//**************************
        daughtermass[0] - sulphur 29 mass
        daughtermass[1] - proton mass

        G4ParticleDefinition*  parent;
        parentmass

//**************************        

00193   //daughters'mass
00194   G4double daughtermass[2]; 
00195   G4double daughtermomentum;
00196   if ( theDaughterMasses )
00197   { 
00198      daughtermass[0]= *(theDaughterMasses);
00199      daughtermass[1] = *(theDaughterMasses+1);
00200   } else {   
00201      daughtermass[0] = daughters[0]->GetPDGMass();
00202      daughtermass[1] = daughters[1]->GetPDGMass();
00203   }
00204   
00205 //  G4double sumofdaughtermass =  daughtermass[0] + daughtermass[1];
00206 
00207   //create parent G4DynamicParticle at rest
00208   G4ParticleMomentum dummy;
00209   G4DynamicParticle * parentparticle = new G4DynamicParticle( parent, dummy, 0.0);
00210 
00211   //create G4Decayproducts  @@GF why dummy parentparticle?
00212   G4DecayProducts *products = new G4DecayProducts(*parentparticle);
00213   delete parentparticle;
00214 
00215   //calculate daughter momentum
00216   daughtermomentum = Pmx(parentmass,daughtermass[0],daughtermass[1]);
00217   G4double costheta = 2.*G4UniformRand()-1.0;
00218   G4double sintheta = std::sqrt((1.0 - costheta)*(1.0 + costheta));
00219   G4double phi  = twopi*G4UniformRand()*rad;
00220   G4ParticleMomentum direction(sintheta*std::cos(phi),sintheta*std::sin(phi),costheta);
00221 
00222   //create daughter G4DynamicParticle
00223   G4double Etotal= std::sqrt(daughtermass[0]*daughtermass[0] + daughtermomentum*daughtermomentum); 
00224   G4DynamicParticle * daughterparticle = new G4DynamicParticle( daughters[0],Etotal, direction*daughtermomentum);
00225   products->PushProducts(daughterparticle);
00226   Etotal= std::sqrt(daughtermass[1]*daughtermass[1] + daughtermomentum*daughtermomentum);
00227   daughterparticle = new G4DynamicParticle( daughters[1],Etotal, direction*(-1.0*daughtermomentum));
00228   products->PushProducts(daughterparticle);
00229 
00230   if (GetVerboseLevel()>1) 
00231     {
00232      G4cout << "G4GeneralPhaseSpaceDecay::TwoBodyDecayIt ";
00233      G4cout << "  create decay products in rest frame " <<G4endl;
00234      products->DumpInfo();
00235     }
00236   return products;

// here I need to convert to lab? I guess 
        rotate with the same angles as costheta etc. from cms to lab
       primGen->AddTrack(pdgType, px, py, pz, vx, vy, vz); //sulpur 
           primGen->AddTrack(pdgType, px, py, pz, vx, vy, vz); //proton
           primGen->AddTrack(pdgType, px, py, pz, vx, vy, vz); //chlorine
00237 }