//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AddVariableAction.h"

class AddPrimarySpeciesAction;

template <>
InputParameters validParams<AddPrimarySpeciesAction>();

class AddPrimarySpeciesAction : public AddVariableAction
{
public:
  AddPrimarySpeciesAction(const InputParameters & params);

  virtual void act() override;

private:
  /// Primary species to add
  const std::vector<NonlinearVariableName> _vars;
  /// Variable scaling
  const std::vector<Real> _scaling;
};