#include "DiffusionCoefficients.h"

DiffusionCoefficient::DiffusionCoefficient() {
	randomSeed = -69;
	boxSize = 256;
	spacing = 0.05 * crpropa::pc;
	brms = 1.0 * crpropa::muG;
	b0 = 0.0 * crpropa::muG;
	nParticles = 1;
	chargeNumber = 1;
	atomicNumber = 1;
	energy = 1e15 * crpropa::eV;
	trajectoryLengthMax = 10 * crpropa::pc;
	maxStep = 0.1 * crpropa::pc;
	doDumpGrid = false;
	nSteps = 0;
	box = NULL;
	prop = NULL;
}

void DiffusionCoefficient::setBrms_muG(double brms) {
	this->brms = brms * crpropa::muG;
}

double DiffusionCoefficient::getBrms_muG() const {
	return (brms / crpropa::muG);
}

void DiffusionCoefficient::setDumpGrid() {
	doDumpGrid = true;
}

void DiffusionCoefficient::unsetDumpGrid() {
	doDumpGrid = false;
}

void DiffusionCoefficient::setSpacing_pc(double spacing) {
	this->spacing = spacing * crpropa::pc;
}

double DiffusionCoefficient::getSpacing_pc() const {
	return (spacing / crpropa::pc);
}

void DiffusionCoefficient::setEnergy_eV(double energy) {
	this->energy = energy * crpropa::eV;
}

double DiffusionCoefficient::getEnergy_eV() const {
	return (energy / crpropa::eV);
}

void DiffusionCoefficient::setB0_muG(double b0) {
	this->b0 = b0 * crpropa::muG;
}

double DiffusionCoefficient::getB0_muG() const {
	return (b0 / crpropa::muG);
}

void DiffusionCoefficient::setTrajectoryLengthMax_pc(double trajectoryLengthMax) {
	this->trajectoryLengthMax = trajectoryLengthMax * crpropa::pc;
}

double DiffusionCoefficient::getTrajectoryLengthMax_pc() const {
	return (trajectoryLengthMax / crpropa::pc);
}

void DiffusionCoefficient::setMaxStep_pc(double maxStep) {
	this->maxStep = maxStep * crpropa::pc;
}

double DiffusionCoefficient::getMaxStep_pc() const {
	return (maxStep / crpropa::pc);
}

void DiffusionCoefficient::setBoxSize(size_t boxSize) {
	this->boxSize = boxSize;
}

size_t DiffusionCoefficient::getBoxSize() const {
	return (boxSize);
}

void DiffusionCoefficient::setNParticles(size_t nParticles) {
	this->nParticles = nParticles;
}

size_t DiffusionCoefficient::getNParticles() const {
	return (nParticles);
}

void DiffusionCoefficient::buildMagneticField(double brms, double b0) {
	setBrms_muG(brms);
	setB0_muG(b0);
	buildMagneticField();
}

void DiffusionCoefficient::loadMagneticField(const std::string& filename) {
	std::cout << "Start loading Magnetic Field..." << std::endl;

	makeMagneticFieldFilename();

	clock_t begin = clock();

	origin.x = 0.0 * crpropa::pc;
	origin.y = 0.0 * crpropa::pc;
	origin.z = 0.0 * crpropa::pc;
	std::cout << "  - Origin: " << origin / crpropa::pc << std::endl;

	int Nx = boxSize;
	int Ny = boxSize;
	int Nz = boxSize;
	int N3 = Nx * Ny * Nz;
	std::cout << "  - Samples: " << Nx << " * " << Ny << " * " << Nz << std::endl;

	size.x = Nx * spacing;
	size.y = Ny * spacing;
	size.z = Nz * spacing;
	std::cout << "  - Spacing: " << spacing / crpropa::pc << " pc " << std::endl;
	std::cout << "  - BoxSize: " << size / crpropa::pc << " pc " << std::endl;

	crpropa::ref_ptr<crpropa::VectorGrid> field = new crpropa::VectorGrid(origin, Nx, Ny, Nz, spacing);

	std::ifstream fin(filename.c_str(), std::ios::binary);
	if (!fin) {
		std::stringstream ss;
		ss << "load MagneticField: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}

	// get length of file and compare to size of grid
	fin.seekg(0, fin.end);
	size_t length = fin.tellg() / sizeof(float);
	fin.seekg(0, fin.beg);

	std::cout << length << "\t" << N3 * 7 << std::endl;

	std::vector<float> density;
	std::vector<float> vx;
	std::vector<float> vy;
	std::vector<float> vz;
	std::vector<float> bx;
	std::vector<float> by;
	std::vector<float> bz;

	float value;

	std::cout << "Reading " << N3 << " characters... " << "\n";
	for (int i = 0; i < N3; ++i) {
		fin.read((char*) &(value), sizeof(float));
		density.push_back(value);
	}
	std::cout << *min_element(density.begin(), density.end()) << "\t" << *max_element(density.begin(), density.end())
			<< "\n";

	std::cout << "Reading " << N3 << " characters... " << "\n";
	for (int i = 0; i < N3; ++i) {
		fin.read((char*) &(value), sizeof(float));
		bx.push_back(value);
	}
	std::cout << "Bx = " << "\t";
	std::cout << bx.size() << "\t";
	std::cout << *min_element(bx.begin(), bx.end()) << "\t" << *max_element(bx.begin(), bx.end()) << "\t";
	std::cout << std::accumulate(bx.begin(), bx.end(), 0.0) / (double) bx.size() << "\n";

	std::cout << "Reading " << N3 << " characters... " << "\n";
	for (int i = 0; i < N3; ++i) {
		fin.read((char*) &(value), sizeof(float));
		by.push_back(value);
	}
	std::cout << "By = " << "\t";
	std::cout << by.size() << "\t";
	std::cout << *min_element(by.begin(), by.end()) << "\t" << *max_element(by.begin(), by.end()) << "\t";
	std::cout << std::accumulate(by.begin(), by.end(), 0.0) / (double) by.size() << "\n";

	std::cout << "Reading " << N3 << " characters... " << "\n";
	for (int i = 0; i < N3; ++i) {
		fin.read((char*) &(value), sizeof(float));
		bz.push_back(value);
	}
	std::cout << "Bz = " << "\t";
	std::cout << bz.size() << "\t";
	std::cout << *min_element(bz.begin(), bz.end()) << "\t" << *max_element(bz.begin(), bz.end()) << "\t";
	std::cout << std::accumulate(bz.begin(), bz.end(), 0.0) / (double) bz.size() << "\n";

	std::cout << "Reading " << N3 << " characters... " << "\n";
	for (int i = 0; i < N3; ++i) {
		fin.read((char*) &(value), sizeof(float));
		vx.push_back(value);
	}
	std::cout << vx.size() << "\t" << *min_element(vx.begin(), vx.end()) << "\t" << *max_element(vx.begin(), vx.end())
			<< "\n";

	std::cout << "Reading " << N3 << " characters... " << "\n";
	for (int i = 0; i < N3; ++i) {
		fin.read((char*) &(value), sizeof(float));
		vy.push_back(value);
	}
	std::cout << vy.size() << "\t" << *min_element(vy.begin(), vy.end()) << "\t" << *max_element(vy.begin(), vy.end())
			<< "\n";

	std::cout << "Reading " << N3 << " characters... " << "\n";
	for (int i = 0; i < N3; ++i) {
		fin.read((char*) &(value), sizeof(float));
		vz.push_back(value);
	}
	std::cout << vz.size() << "\t" << *min_element(vz.begin(), vz.end()) << "\t" << *max_element(vz.begin(), vz.end())
			<< "\n";

	if (fin)
		std::cout << "all characters read successfully." << "\n";
	else
		std::cout << "error: only " << fin.gcount() << " could be read";

	int counter = 0;
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++) {
				crpropa::Vector3f &b = field->get(ix, iy, iz);
				b.x = bx[counter] * crpropa::gauss;
				b.y = by[counter] * crpropa::gauss;
				b.z = bz[counter] * crpropa::gauss;
				counter++;
			}
	fin.close();

	magneticField = new crpropa::MagneticFieldGrid(field);

	clock_t end = clock();

	std::cout << "  - Elapsed secs = " << double(end - begin) / CLOCKS_PER_SEC << std::endl;
}

void DiffusionCoefficient::buildMagneticField() {
	std::cout << "Start building Magnetic Field..." << std::endl;

	makeMagneticFieldFilename();

	clock_t begin = clock();

	crpropa::Random::seedThreads(randomSeed);

	origin.x = 0.0 * crpropa::pc;
	origin.y = 0.0 * crpropa::pc;
	origin.z = 0.0 * crpropa::pc;
	std::cout << "  - Origin: " << origin / crpropa::pc << std::endl;

	int Nx = boxSize;
	int Ny = boxSize;
	int Nz = boxSize;
	std::cout << "  - Samples: " << Nx << " * " << Ny << " * " << Nz << std::endl;

	size.x = Nx * spacing;
	size.y = Ny * spacing;
	size.z = Nz * spacing;
	std::cout << "  - Spacing: " << spacing / crpropa::pc << " pc " << std::endl;
	std::cout << "  - BoxSize: " << size / crpropa::pc << " pc " << std::endl;

	crpropa::ref_ptr<crpropa::VectorGrid> field = new crpropa::VectorGrid(origin, Nx, Ny, Nz, spacing);

	std::cout << "  - Creating random turbulent field" << std::endl;
	std::cout << "  - Brms: " << brms / crpropa::muG << " muG" << std::endl;
	std::cout << "  - Breg: " << b0 / crpropa::muG << " muG" << std::endl;

	double kMin = 2.0 / (double) Nx;
	double kMax = 0.5;
	double lMin = spacing / kMax;
	double lMax = spacing / kMin;

	std::cout << "  - Turbulent range: " << lMin / crpropa::pc << " - " << lMax / crpropa::pc << " pc" << std::endl;

	double alpha = -11.0 / 3.0;

	std::cout << "  - Spectral index, <B^2(k)> ~ k^n, n:  " << alpha << std::endl;

	crpropa::initTurbulence(field, brms, lMin, lMax, alpha, 42, false, 0);

	if (doDumpGrid) {
		crpropa::dumpGrid(field, magneticFieldFilename);
	}

	magneticField = new crpropa::MagneticFieldGrid(field);

	clock_t end = clock();

	std::cout << "  - Elapsed secs = " << double(end - begin) / CLOCKS_PER_SEC << std::endl;
}

void DiffusionCoefficient::buildCandidates(size_t nParticles, double energy) {
	setEnergy_eV(energy);
	setNParticles(nParticles);
	buildCandidates();
}

void DiffusionCoefficient::buildCandidates() {
	std::cout << "Start building " << nParticles << " Candidates..." << std::endl;

	for (size_t i = 0; i < nParticles; ++i) {
		crpropa::Random & random = crpropa::Random::instance();
		crpropa::Vector3d position(random.rand(), random.rand(), random.rand());
		position = position * size + origin;
		crpropa::Vector3d direction = crpropa::Vector3d(random.randUniform(-1, 1), random.randUniform(-1, 1),
				random.randUniform(-1, 1));
		direction = direction.getUnitVector();
		int id = crpropa::nucleusId(atomicNumber, chargeNumber);

		crpropa::ref_ptr<crpropa::Candidate> candidate = new crpropa::Candidate(id, energy, position, direction);
		//std::cout << candidate->getDescription() << "\n";

		particles.push_back(candidate);
	}
}

void DiffusionCoefficient::buildAlgorithm() {
	std::cout << "Start building Propagation..." << std::endl;

	double epsilon = 1e-4; // MinStep_pc
	std::cout << "  - Epsilon: " << epsilon << std::endl;

	double minStep = epsilon * crpropa::pc; // MinStep_pc
	std::cout << "  - Minimum step: " << minStep / crpropa::pc << " pc" << std::endl;

	//double maxStep = 1e-1 * crpropa::pc; // MaxStep_pc
	std::cout << "  - Maximum step: " << maxStep / crpropa::pc << " pc" << std::endl;

	if (minStep >= maxStep)
		throw std::runtime_error(" --> Maximum step must be larger than minimum step");

	if ((size.x == 0) or (size.y == 0) or (size.z == 0))
		throw std::runtime_error(" --> Environment boundaries not set");

	std::cout << "  - Lower bounds: " << origin / crpropa::pc << " pc" << std::endl;
	std::cout << "  - Upper bounds: " << (origin + size) / crpropa::pc << " pc" << std::endl;

	box = new crpropa::PeriodicBox(origin, size);

	std::cout << box->getDescription() << std::endl;

	prop = new crpropa::PropagationCK(magneticField, epsilon, minStep, maxStep);

	std::cout << prop->getDescription() << std::endl;
}

void DiffusionCoefficient::makeMagneticFieldFilename() {
	std::stringstream ss;
	ss << "magneticField_" << brms / crpropa::muG << "_" << b0 / crpropa::muG << "_" << boxSize << ".bin";
	magneticFieldFilename = ss.str();
}

void DiffusionCoefficient::makeTrajectoryFilename() {
	std::stringstream ss;
	ss << "trajectory_" << nParticles << "_" << energy / crpropa::PeV << "_" << brms / crpropa::muG << "_"
			<< b0 / crpropa::muG << "_" << boxSize << ".txt";
	trajectoryFilename = ss.str();
}

void DiffusionCoefficient::buildTrajectories() {
	for (size_t i = 0; i < nParticles; ++i) {
		Trajectory* t = new Trajectory(i);
		trajs.push_back(t);
	}
}

void DiffusionCoefficient::run(double trajectoryLengthMax) {
	setTrajectoryLengthMax_pc(trajectoryLengthMax);
	run();
}
;

void DiffusionCoefficient::run() {
	std::cout << "Run..." << std::endl;
	std::cout << "  - trajectoryLengthMax = " << trajectoryLengthMax / crpropa::pc << " pc" << std::endl;

	clock_t begin = clock();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(OMP_NUM_THREADS)
#endif
	for (size_t i = 0; i < particles.size(); ++i) {
		size_t counter = 0;
		while (particles[i]->getTrajectoryLength() < trajectoryLengthMax) {
			prop->process(particles[i]);
			box->process(particles[i]);
			if (counter % 2 == 0) {
				trajs[i]->addTrajectoryLength(particles[i]->getTrajectoryLength());
				trajs[i]->addPosition(particles[i]->getCurrentPosition(), particles[i]->getCreatedPosition());
				trajs[i]->addCorrelations(particles[i]->getCurrentPosition(), particles[i]->getCreatedPosition(),
						particles[i]->getCurrentDirection());
			}
			counter++;
		}
	}

	nSteps = trajs[0]->getLenghtVector();

	clock_t end = clock();
	std::cout << "  - Elapsed secs = " << double(end - begin) / CLOCKS_PER_SEC << std::endl;
	std::cout << "  - Trajectory Vector size = " << nSteps << std::endl;
}

void DiffusionCoefficient::calcOutput() {
	for (size_t n = 0; n < nSteps; ++n) {
		double cdt = trajs[0]->getCT(n);
		crpropa::Vector3d dist;
		crpropa::Vector3d dist2;
		crpropa::Vector3d r;
		for (size_t i = 0; i < nParticles; ++i) {
			dist += trajs[i]->getPosition(n);
			dist2 += trajs[i]->getPosition2(n);
			r += trajs[i]->getR(n);
		}
		outputParams o = { cdt, dist / (double) nParticles, dist2 / (double) nParticles, r / (double) nParticles };
		output.push_back(o);
	}
}

void DiffusionCoefficient::writeOutput() {
	std::cout << "Writing Output..." << std::endl;

	makeTrajectoryFilename();
	calcOutput();

	std::ofstream trajectoryFile(trajectoryFilename.c_str());
	trajectoryFile << std::scientific << std::setprecision(5);
	for (size_t i = 0; i < output.size(); ++i) {
		trajectoryFile << output[i].cdt / crpropa::pc << "\t";
		trajectoryFile << output[i].distAverage.getX() / crpropa::pc << "\t";
		trajectoryFile << output[i].distAverage.getY() / crpropa::pc << "\t";
		trajectoryFile << output[i].distAverage.getZ() / crpropa::pc << "\t";
		trajectoryFile << output[i].dist2Average.getX() / pow(crpropa::pc, 2.0) << "\t";
		trajectoryFile << output[i].dist2Average.getY() / pow(crpropa::pc, 2.0) << "\t";
		trajectoryFile << output[i].dist2Average.getZ() / pow(crpropa::pc, 2.0) << "\t";
		trajectoryFile << output[i].RijAverage.getX() / crpropa::pc << "\t";
		trajectoryFile << output[i].RijAverage.getY() / crpropa::pc << "\t";
		trajectoryFile << output[i].RijAverage.getZ() / crpropa::pc << "\n";
	}
	trajectoryFile.close();
}

DiffusionCoefficient::~DiffusionCoefficient() {
	std::cerr << "Deleting DiffusionCoefficient:" << std::endl;

	if (output.size() > 0) {
		std::cerr << "  - clearing vector<outputParams> output" << std::endl;
		output.clear();
	}
	if (particles.size() > 0) {
		std::cerr << "  - clearing vector<ref_ptr<Candidate> > particles" << std::endl;
		particles.clear();
	}
	if (trajs.size() > 0) {
		std::cerr << "  - clearing vector<ref_ptr<Trajectory> > trajs" << std::endl;
		trajs.clear();
	}
	if (box) {
		std::cerr << "  - deleting PeriodicBox*" << std::endl;
		delete box;
	}
	if (prop) {
		std::cerr << "  - deleting PropagationCK*" << std::endl;
		delete prop;
	}
}
