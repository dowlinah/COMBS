//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// STL includes
#include <fstream>

// MOOSE includes
#include "CSVReader.h"
#include "MooseUtils.h"

registerMooseObject("MooseApp", CSVReader);

defineLegacyParams(CSVReader);

InputParameters
CSVReader::validParams()
{
  InputParameters params = GeneralVectorPostprocessor::validParams();
  params.addClassDescription(
      "Converts columns of a CSV file into vectors of a VectorPostprocessor.");
  params.addRequiredParam<FileName>("csv_file",
                                    "The name of the CSV file to read. Currently, with "
                                    "the exception of the header row, only numeric "
                                    "values are supported.");
  params.addParam<bool>("header",
                        "When true it is assumed that the first row contains column headers, these "
                        "headers are used as the VectorPostprocessor vector names. If false the "
                        "file is assumed to contain only numbers and the vectors are named "
                        "automatically based on the column number (e.g., 'column_0000', "
                        "'column_0001'). If not supplied the reader attempts to auto detect the "
                        "headers.");
  params.addParam<std::string>("delimiter",
                               "The column delimiter. Despite the name this can read files "
                               "separated by delimiter other than a comma. If this options is "
                               "omitted it will read comma or space separated files.");
  params.addParam<bool>(
      "ignore_empty_lines", true, "When true new empty lines in the file are ignored.");
  params.set<ExecFlagEnum>("execute_on", true) = EXEC_INITIAL;

  // The value from this VPP is naturally already on every processor
  // TODO: Make this not the case!  See #11415
  params.set<bool>("_is_broadcast") = true;

  return params;
}

CSVReader::CSVReader(const InputParameters & params)
  : GeneralVectorPostprocessor(params), _csv_reader(getParam<FileName>("csv_file"), &_communicator)
{
  _csv_reader.setIgnoreEmptyLines(getParam<bool>("ignore_empty_lines"));
  if (isParamValid("header"))
    _csv_reader.setHeaderFlag(getParam<bool>("header")
                                  ? MooseUtils::DelimitedFileReader::HeaderFlag::ON
                                  : MooseUtils::DelimitedFileReader::HeaderFlag::OFF);
  if (isParamValid("delimiter"))
    _csv_reader.setDelimiter(getParam<std::string>("delimiter"));
}

void
CSVReader::initialize()
{
  _csv_reader.read();
  for (auto & name : _csv_reader.getNames())
    if (_column_data.find(name) == _column_data.end())
      _column_data[name] = &declareVector(name);
}

void
CSVReader::execute()
{
  const std::vector<std::string> & names = _csv_reader.getNames();
  const std::vector<std::vector<double>> & data = _csv_reader.getData();
  for (std::size_t i = 0; i < _column_data.size(); ++i)
    _column_data[names[i]]->assign(data[i].begin(), data[i].end());
}
