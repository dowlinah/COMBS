//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FluxBC.h"

class HeatConductionBC;

template <>
InputParameters validParams<HeatConductionBC>();

/**
 *
 */
class HeatConductionBC : public FluxBC
{
public:
  static InputParameters validParams();

  HeatConductionBC(const InputParameters & parameters);
  virtual ~HeatConductionBC();

protected:
  virtual RealGradient computeQpFluxResidual();
  virtual RealGradient computeQpFluxJacobian();

  const MaterialProperty<Real> & _k;
};

