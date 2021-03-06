//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DarcyPressure.h"

registerADMooseObject("DarcyThermoMechApp", DarcyPressure);

template <ComputeStage compute_stage>
InputParameters
DarcyPressure<compute_stage>::validParams()
{
  InputParameters params = ADKernel<compute_stage>::validParams();
  params.addClassDescription("Compute the diffusion term for Darcy pressure ($p$) equation: "
                             "$-\\nabla \\cdot \\frac{\\mathbf{K}}{\\mu} \\nabla p = 0$");

  // Add a required parameter.  If this isn't provided in the input file MOOSE will error.
  params.addRequiredParam<Real>("permeability", "The permeability ($\\mathrm{K}$) of the fluid.");

  // Add a parameter with a default value; this value can be overridden in the input file.
  params.addParam<Real>(
      "viscosity",
      7.98e-4,
      "The viscosity ($\\mu$) of the fluid in Pa, the default is for water at 30 degrees C.");
  return params;
}

template <ComputeStage compute_stage>
DarcyPressure<compute_stage>::DarcyPressure(const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters),

    // Get the parameters from the input file
    _permeability(getParam<Real>("permeability")),
    _viscosity(getParam<Real>("viscosity"))
{
}

template <ComputeStage compute_stage>
ADReal
DarcyPressure<compute_stage>::computeQpResidual()
{
  return (_permeability / _viscosity) * _grad_test[_i][_qp] * _grad_u[_qp];
}
