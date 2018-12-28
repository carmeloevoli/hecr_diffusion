#include "Trajectory.h"

Trajectory::Trajectory(size_t id) {
	this->id = id;
}

void Trajectory::addPosition(crpropa::Vector3d currentPosition,
		crpropa::Vector3d createdPosition) {
	positions.push_back(currentPosition - createdPosition);
}

void Trajectory::addTrajectoryLength(double trajectoryLength) {
	cdt.push_back(trajectoryLength);
}

void Trajectory::addCorrelations(crpropa::Vector3d currentPosition,
		crpropa::Vector3d createdPosition, crpropa::Vector3d currentDirection) {
	double Rxx = currentDirection.getX()
			* (currentPosition.getX() - createdPosition.getX());
	double Rzz = currentDirection.getZ()
			* (currentPosition.getZ() - createdPosition.getZ());
	double Ryz = currentDirection.getY()
			* (currentPosition.getZ() - createdPosition.getZ());

	crpropa::Vector3d R(Rxx, Rzz, Ryz);

	Rij.push_back(R);
}

size_t Trajectory::getLenghtVector() const {
	return (positions.size());
}

int Trajectory::getId() const {
	return (id);
}

crpropa::Vector3d Trajectory::getR(size_t i) const {
	return (Rij[i]);
}

crpropa::Vector3d Trajectory::getPosition(size_t i) const {
	return (positions[i]);
}

crpropa::Vector3d Trajectory::getPosition2(size_t i) const {
	return (positions[i] * positions[i]);
}

double Trajectory::getCT(size_t i) const {
	return (cdt[i]);
}
