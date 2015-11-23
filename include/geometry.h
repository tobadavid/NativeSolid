#pragma once

#include "config.h"
#include <vector>

class Point{

protected:
  double* Coord;
  double* Coord_n;

public:
  Point();
  ~Point();
  double* GetCoord() const;
  double* GetCoord_n() const;
  void PrintCoord() const;
  void SetCoord(double* newCoord);
  void SetCoord_n(double* newCoord);
  void UpdateCoord();
};

/*class Vertex{

protected:
  double* Coord;

public:
  Vertex();
  ~Vertex();
  double* GetCoord() const;
  void SetCoord(double* newCoord);

}*/

/*class Element{

protected:
  Point** pointElem;
  unsigned short nPointElem;

public:
  Element();
  ~Element();
}*/

/*class Marker{

protected:
  std::string MarkerTag;
  unsigned long nElem;

public:
  Marker(std::string Tag, unsigned long nElem_value);
  ~Marker();
  Vertex** vertex;
  Element** element;

}*/

class Geometry{

protected:
  unsigned short nDim;
  unsigned long nElem, nPoint, nMarkers;
  unsigned long* nElemMarker;

public:
  Geometry(Config* config);
  ~Geometry();
  unsigned long GetnMarkers() const;
  bool GetMarkersMoving(unsigned long iMarker) const;
  void UpdateGeometry();
  unsigned long** vertex;  //vertex[iMarker][iPoint]
  unsigned long* nVertex;  //nVertex[iMarker]
  bool* markersMoving;
  Point** node;         //node[iPöint] returns a pointeur to a Point object so e.g. node[iPoint]->PrintCoord();
};

bool isInVec(std::vector<int> const& inputVector, int dummyInt);
void vecCopy(std::vector<int> const& sourceVector, unsigned long* destinationTab);
