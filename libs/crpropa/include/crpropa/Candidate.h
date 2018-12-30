#ifndef CRPROPA_CANDIDATE_H
#define CRPROPA_CANDIDATE_H

#include "crpropa/ParticleState.h"
#include "crpropa/Referenced.h"
#include "crpropa/AssocVector.h"
#include "crpropa/Variant.h"

#include <vector>
#include <map>
#include <sstream>
#include <stdint.h>

namespace crpropa {
/**
 * \addtogroup Core
 * @{
 */

/**
 @class Candidate
 @brief All information about the cosmic ray.

 The Candidate is a passive object, that holds the information about the state
 of the cosmic ray and the simulation itself.
 */
class Candidate: public Referenced {
public:
	ParticleState source; /**< Particle state at the source */
	ParticleState created; /**< Particle state of parent particle at the time of creation */
	ParticleState current; /**< Current particle state */
	ParticleState previous; /**< Particle state at the end of the previous step */

	typedef Loki::AssocVector<std::string, Variant> PropertyMap;
	PropertyMap properties; /**< Map of property names and their values. */

private:
	bool active; /**< Active status */
	double trajectoryLength; /**< Comoving distance [m] the candidate has traveled so far */
	double currentStep; /**< Size of the currently performed step in [m] comoving units */
	double nextStep; /**< Proposed size of the next propagation step in [m] comoving units */

	static uint64_t nextSerialNumber;
	uint64_t serialNumber;

public:
	Candidate(
		int id = 0,
		double energy = 0,
		Vector3d position = Vector3d(0, 0, 0),
		Vector3d direction = Vector3d(-1, 0, 0)
	);

	/**
	 Creates a candidate, initializing the Candidate::source, Candidate::created,
	 Candidate::previous and Candidate::current state with the argument.
	 */
	Candidate(const ParticleState &state);

	const Vector3d &getCurrentPosition() const {
		return current.getPosition();
	}

	const Vector3d &getCreatedPosition() const {
		return created.getPosition();
	}

	const Vector3d &getCurrentDirection() const {
		return current.getDirection();
	}

	bool isActive() const;
	void setActive(bool b);

	void setTrajectoryLength(double length);
	double getTrajectoryLength() const;

	/**
	 Sets the current step and increases the trajectory length accordingly.
	 Only the propagation module should use this.
	 */
	void setCurrentStep(double step);
	double getCurrentStep() const;

	/**
	 Sets the proposed next step.
	 Only the propagation module should use this.
	 */
	void setNextStep(double step);
	double getNextStep() const;

	/**
	 Make a bid for the next step size: the lowest wins.
	 */
	void limitNextStep(double step);

	void setProperty(const std::string &name, const Variant &value);
	const Variant &getProperty(const std::string &name) const;
	bool removeProperty(const std::string &name);
	bool hasProperty(const std::string &name) const;

	std::string getDescription() const;

	/** Unique (inside process) serial number (id) of candidate */
	uint64_t getSerialNumber() const;
	void setSerialNumber(const uint64_t snr);

	/** Serial number of candidate at source*/
	uint64_t getSourceSerialNumber() const;

	/** Serial number of candidate at creation */
	uint64_t getCreatedSerialNumber() const;

	/** Set the next serial number to use */
	static void setNextSerialNumber(uint64_t snr);

	/** Get the next serial number that will be assigned */
	static uint64_t getNextSerialNumber();

	/**
	 Create an exact clone of candidate
	 @param recursive	recursively clone and add the secondaries
	 */
	ref_ptr<Candidate> clone(bool recursive = false) const;

	/**
	 Copy the source particle state to the current state
	 and activate it if inactive, e.g. restart it
	*/
	void restart();
};

/** @}*/
} // namespace crpropa

#endif // CRPROPA_CANDIDATE_H
