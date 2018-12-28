#ifndef __TRAJECTORY_H
#define __TRAJECTORY_H

#include <vector>

#include "crpropa/Referenced.h"
#include "crpropa/Vector3.h"

class Trajectory: public crpropa::Referenced {
private:
	int id;
	std::vector<double> cdt;
	std::vector<crpropa::Vector3d> positions;
	std::vector<crpropa::Vector3d> Rij;

public:
	Trajectory(size_t id);
	~Trajectory() {
	}

	void addPosition(crpropa::Vector3d currentPosition,
			crpropa::Vector3d createdPosition);
	void addTrajectoryLength(double trajectoryLength);
	void addCorrelations(crpropa::Vector3d currentPosition,
			crpropa::Vector3d createdPosition,
			crpropa::Vector3d currentDirection);
	size_t getLenghtVector() const;
	int getId() const;
	double getCT(size_t i) const;
	crpropa::Vector3d getPosition(size_t i) const;
	crpropa::Vector3d getPosition2(size_t i) const;
	crpropa::Vector3d getR(size_t i) const;
};

#endif // __TRAJECTORY_
