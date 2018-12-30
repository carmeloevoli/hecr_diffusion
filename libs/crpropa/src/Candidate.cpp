#include "crpropa/Candidate.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Units.h"

#include <stdexcept>

namespace crpropa {

Candidate::Candidate(int id, double E, Vector3d pos, Vector3d dir) :
		trajectoryLength(0), currentStep(0), nextStep(0), active(true) {
	ParticleState state(id, E, pos, dir);
	source = state;
	created = state;
	previous = state;
	current = state;

#if defined(OPENMP_3_1)
#pragma omp atomic capture
	{	serialNumber = nextSerialNumber++;}
#elif defined(__GNUC__)
	{
		serialNumber = __sync_add_and_fetch(&nextSerialNumber, 1);
	}
#else
#pragma omp critical
	{	serialNumber = nextSerialNumber++;}
#endif

}

Candidate::Candidate(const ParticleState &state) :
		source(state), created(state), current(state), previous(state), trajectoryLength(0), currentStep(0), nextStep(
				0), active(true) {

#if defined(OPENMP_3_1)
#pragma omp atomic capture
	{	serialNumber = nextSerialNumber++;}
#elif defined(__GNUC__)
	{
		serialNumber = __sync_add_and_fetch(&nextSerialNumber, 1);
	}
#else
#pragma omp critical
	{	serialNumber = nextSerialNumber++;}
#endif

}

bool Candidate::isActive() const {
	return active;
}

void Candidate::setActive(bool b) {
	active = b;
}

double Candidate::getTrajectoryLength() const {
	return trajectoryLength;
}

double Candidate::getCurrentStep() const {
	return currentStep;
}

double Candidate::getNextStep() const {
	return nextStep;
}

void Candidate::setTrajectoryLength(double a) {
	trajectoryLength = a;
}

void Candidate::setCurrentStep(double lstep) {
	currentStep = lstep;
	trajectoryLength += lstep;
}

void Candidate::setNextStep(double step) {
	nextStep = step;
}

void Candidate::limitNextStep(double step) {
	nextStep = std::min(nextStep, step);
}

void Candidate::setProperty(const std::string &name, const Variant &value) {
	properties[name] = value;
}

const Variant &Candidate::getProperty(const std::string &name) const {
	PropertyMap::const_iterator i = properties.find(name);
	if (i == properties.end())
		throw std::runtime_error("Unknown candidate property: " + name);
	return i->second;
}

bool Candidate::removeProperty(const std::string& name) {
	PropertyMap::iterator i = properties.find(name);
	if (i == properties.end())
		return false;
	properties.erase(i);
	return true;
}

bool Candidate::hasProperty(const std::string &name) const {
	PropertyMap::const_iterator i = properties.find(name);
	if (i == properties.end())
		return false;
	return true;
}


std::string Candidate::getDescription() const {
	std::stringstream ss;
	ss << "CosmicRay source: " << source.getDescription() << "\n";
	ss << "     ... current: " << current.getDescription();
	return ss.str();
}

ref_ptr<Candidate> Candidate::clone(bool recursive) const {
	ref_ptr<Candidate> cloned = new Candidate;
	cloned->source = source;
	cloned->created = created;
	cloned->current = current;
	cloned->previous = previous;

	cloned->properties = properties;
	cloned->active = active;
	cloned->trajectoryLength = trajectoryLength;
	cloned->currentStep = currentStep;
	cloned->nextStep = nextStep;
	return cloned;
}

uint64_t Candidate::getSerialNumber() const {
	return serialNumber;
}

void Candidate::setSerialNumber(const uint64_t snr) {
	serialNumber = snr;
}

uint64_t Candidate::getSourceSerialNumber() const {
		return serialNumber;
}

uint64_t Candidate::getCreatedSerialNumber() const {
		return serialNumber;
}

void Candidate::setNextSerialNumber(uint64_t snr) {
	nextSerialNumber = snr;
}

uint64_t Candidate::getNextSerialNumber() {
	return nextSerialNumber;
}

uint64_t Candidate::nextSerialNumber = 0;

void Candidate::restart() {
	setActive(true);
	setTrajectoryLength(0);
	previous = source;
	current = source;
}

} // namespace crpropa
