#include <cmath>
#include <iostream>

#include "DiffusionCoefficients.h"

#define NGRID 128
#define ENERGY 1e15
#define NPARTICLES 10
#define BOXSIZE 120

int main(int argc, char** argv) {

	try {

		DiffusionCoefficient D;
		D.setBoxSize(NGRID);
		D.setSpacing_pc((double) BOXSIZE / (double) NGRID);
		D.setBrms_muG(1.0);
		D.setB0_muG(1.0);
		D.setDumpGrid();
		D.setMaxStep_pc(0.01 * (ENERGY / 1e15));
		D.setTrajectoryLengthMax_pc(1e4);
		D.setNParticles(NPARTICLES);
		D.setEnergy_eV(ENERGY);
		D.buildMagneticField();
		D.buildCandidates();
		D.buildTrajectories();
		D.buildAlgorithm();
		D.run();
		D.writeOutput();
	} // try

	catch (std::exception &e) {
		std::cerr << "Exception: " << e.what() << std::endl;
		return (1);
	}

	return (0);
}

//int main(int argc, char** argv) {

//  Brandom("box/3D_fields_10_64_1_31.bin",64,10.0);

//  if (argc < 6) {
//  cout<<"Usage: ./diffusion NBFieldsMax NParticlesMax gamma turbulenceLevel turbulenceDimensions <filename>"<<endl;
//  cout<<"Example: ./a.out 1 1 1e5 2 3 test"<<endl;
//  return -1;
//}

//     const double rL = A / Z * PROTON_MASS * sqrt(gamma * gamma - 1.0) * TMath::C()/(electric_charge*1e-10)/parsec; 
//  const int NBFieldsMax = atoi(argv[1]); // Max number of field configurations
//   const int NParticlesMax = atoi(argv[2]); // Max number of particles
//   const double gamma = atof(argv[3]); // EnergyProton/MassProtons. It's better to set this first.
//   const double turbulenceLevel = atof(argv[4]); // sigma2
//   const int turbulenceDimensions = atoi(argv[5]);

//   // Maximal correlation lenght (100 pc) expressed in units of r_L

//   const double rL =  A/Z*PROTON_MASS*sqrt(gamma*gamma-1.0)*TMath::C()/(electric_charge*1e-10)/parsec; 
//   const double Lmax = 100.0/rL;
//   const double LcrL = Lc/rL;
//   const long NLarmorTimes = /*16000.0*/1000.0;
//   const double tEnd = (double)NLarmorTimes*(rL/Clight); // final time in sec
//   const double dt = 0.5*Lmin*(rL/Clight);
//   const int NSteps = int(tEnd/dt)+1;

//   vector< vector<double> > results(NSteps, vector<double>(7,0));

//   cout<<"... dt = "<<dt<<endl;
//   cout<<"... tEnd = "<<tEnd<<endl;
//   cout<<"... rL = "<<rL<<endl;
//   cout<<"... Lc*rL = "<<LcrL<<endl;

//   // Time counters...
//   TStopwatch time_s; 
//   time_s.Start();

//   TNtupleD* nt = new TNtupleD("nt", "Particle correlation data", "t:xx0:yy0:zz0:vxx:vxy:vyy:vzz");

// #ifdef _OPENMP
// #pragma omp parallel //for schedule(dynamic) num_threads(OMP_NUM_THREADS)
// #endif
//   {
//     for ( int i=0;i<NBFieldsMax;++i ){ // Loop over magnetic field configurations

//       Brandom Br();//KOLMOGOROV, Lmin, Lmax, LcrL, turbulenceLevel, turbulenceDimensions, 0);

//       TVector3 x(1./3.,1./3.,1./3.);
//       //Br.Set(x);
//       //cout<<Br.Bx()<<"\t"<<Br.By()<<"\t"<<Br.Bz()<<endl;

// #ifdef _OPENMP
// #pragma omp for schedule(dynamic) private(Br) nowait
// #endif
//       for ( int l=0;l<NParticlesMax;++l ){ // Loop over the particles    

// #ifdef _OPENMP
// #pragma omp critical
// #endif
// 	cout<<"iField = "<<i<<"  iParticle = "<<l<<endl;
// 	int id = i*NParticlesMax+l;
// 	Particle part(0, gamma, tEnd, dt, rL, Br, id);
// 	vector< vector<double> > output(NSteps, vector<double>(7,0));

// 	part.Evolve(output);

// #ifdef _OPENMP
// #pragma omp critical
// #endif
// 	{
// 	  for ( int itime=0;itime<NSteps;itime++ ){
// 	    for ( int iout=0;iout<7;iout++ ) 
// 	      results[itime][iout] += output[itime][iout]/double(NBFieldsMax*NParticlesMax);
// 	  } 
// 	}

//       }
//     }
//   }

// #ifdef ROOTOUTPUT
//   TFile* rootfile = TFile::Open(argv[argc-1], "RECREATE");

//   TNtupleD* header = new TNtupleD("header", "Header", "rL:energy");
//   header->Fill(rL, PROTON_MASS_GEV*gamma);
//   header->Write();

//   for ( int itime=0;itime<NSteps;itime++ ) 
//     nt->Fill(double(itime)*dt, 
// 	     results[itime][0],
// 	     results[itime][1],
// 	     results[itime][2],
// 	     results[itime][3],
// 	     results[itime][4],
// 	     results[itime][5],
// 	     results[itime][6]
// 	     );

//   nt->Write();
//   rootfile->Close();
// #endif

// #ifdef ASCIIOUTPUT
//   ofstream asciifile("test.txt");

//   for ( int itime=0;itime<NSteps;itime++ )
//     asciifile<<double(itime)*dt<<"\t"
// 	     <<results[itime][0]<<"\t"
// 	     <<results[itime][1]<<"\t"
// 	     <<results[itime][2]<<"\t"
// 	     <<results[itime][3]<<"\t"
// 	     <<results[itime][4]<<"\t"
// 	     <<results[itime][5]<<"\t"
// 	     <<results[itime][6]<<endl;

//   asciifile.close();
// #endif  

//   time_s.Stop();

//   cout<<"... Ended in "<<time_s.RealTime()<<" s."<<endl;
//   cout<<"... Actual CPU time "<<time_s.CpuTime()<<" s."<<endl;
//   cout<<"... Time per particle "<<time_s.RealTime()/double(NBFieldsMax*NParticlesMax)<<" s."<<endl;

//  return 0;
//}
