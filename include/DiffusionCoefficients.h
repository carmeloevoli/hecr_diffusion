#ifndef __DIFFUSION_COEFFICIENT_H
#define __DIFFUSION_COEFFICIENT_H

#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include "crpropa/Candidate.h"
#include "crpropa/Grid.h"
#include "crpropa/GridTools.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Random.h"
#include "crpropa/Vector3.h"
#include "crpropa/Units.h"
#include "crpropa/magneticField/MagneticFieldGrid.h"
#include "crpropa/module/Boundary.h"
#include "crpropa/module/PropagationCK.h"

#include "Trajectory.h"

#ifdef _OPENMP
#include <omp.h>
#define OMP_NUM_THREADS 16
#endif

struct outputParams {
	double cdt;
	crpropa::Vector3d distAverage;
	crpropa::Vector3d dist2Average;
	crpropa::Vector3d RijAverage;
};

class DiffusionCoefficient {
private:

	bool doDumpGrid;

	size_t nParticles;
	size_t boxSize;
	size_t nSteps;
	int randomSeed;
	int chargeNumber;
	int atomicNumber;
	double energy;
	double spacing;
	double brms;
	double b0;
	double trajectoryLengthMax;
	double maxStep;

	std::string magneticFieldFilename;
	std::string trajectoryFilename;

	std::vector<crpropa::ref_ptr<crpropa::Candidate> > particles;
	std::vector<crpropa::ref_ptr<Trajectory> > trajs;
	std::vector<outputParams> output;

	crpropa::Vector3d origin;
	crpropa::Vector3d size;
	crpropa::ref_ptr<crpropa::MagneticField> magneticField;
	crpropa::PeriodicBox* box;
	crpropa::PropagationCK* prop;

	void makeMagneticFieldFilename();
	void makeTrajectoryFilename();

public:

	DiffusionCoefficient();
	~DiffusionCoefficient();

	void setDumpGrid();
	void unsetDumpGrid();

	void setBrms_muG(double brms);
	double getBrms_muG() const;

	void setB0_muG(double b0);
	double getB0_muG() const;

	void setSpacing_pc(double spacing);
	double getSpacing_pc() const;

	void setEnergy_eV(double energy);
	double getEnergy_eV() const;

	void setTrajectoryLengthMax_pc(double TrajectoryLengthMax);
	double getTrajectoryLengthMax_pc() const;

	void setMaxStep_pc(double maxStep);
	double getMaxStep_pc() const;

	void setBoxSize(size_t boxSize);
	size_t getBoxSize() const;

	void setNParticles(size_t nParticles);
	size_t getNParticles() const;

	void buildMagneticField(double brms, double b0);
	void buildMagneticField();

	void buildCandidates(size_t nParticles, double energy);
	void buildCandidates();
	void loadMagneticField(const std::string& filename);

	void buildAlgorithm();

	void buildTrajectories();

	void run(double trajectoryLengthMax);
	void run();

	void calcOutput();
	void writeOutput();
};

#endif
