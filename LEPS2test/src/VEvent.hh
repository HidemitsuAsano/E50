/*
  VEvent.hh

  2012/5  K.Shirotori
*/

#ifndef VEVENT_H
#define VEVENT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

class VEvent
{

public:
  VEvent();
  virtual ~VEvent() = 0;

  virtual bool ProcessingBegin() = 0;
  virtual bool ProcessingEnd() = 0;
  virtual bool ProcessingNormal( std::ifstream & ) = 0;

};

#endif
