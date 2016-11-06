/*
  VisManager.hh

  2012/5  K.Shirotori
*/

#ifndef VIS_MANAGER_H
#define VIS_MANAGER_H 1

#include "G4VisManager.hh"

class VisManager : public G4VisManager {
public:
  VisManager();
  ~VisManager();

private:
  void RegisterGraphicsSystems();
  void RegisterModelFactories();
};

#endif
