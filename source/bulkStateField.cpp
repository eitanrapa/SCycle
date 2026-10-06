#include "bulkStateField.hpp"
#include "powerLaw.hpp"
#include <fstream>
#include <cmath>
#include <algorithm>

#define FILENAME "bulkStateField.cpp"

using namespace std;


BulkStateField::BulkStateField(Domain& D, const string& name, const string& prefix, const string& symbol)
: _D(&D), _file(D._file), _delim(D._delim), _inputDir(D._inputDir), _outputDir(D._outputDir),
  _name(name), _prefix(prefix), _symbol(symbol), _type("off"),
  _state(NULL), _rate(NULL), _e(NULL),
  _eTest(-1.0), _eTestAmp(0.0), _eTestPeriod(1.0)
{}


BulkStateField::~BulkStateField()
{
  VecDestroy(&_state);
  VecDestroy(&_rate);
  VecDestroy(&_e);
}


void BulkStateField::readLines(vector<string>& vars, vector<string>& rhss, vector<string>& rhsFulls) const
{
  ifstream infile(_file);
  string line;
  while (getline(infile, line)) {
    size_t pos = line.find(_delim);
    string var = line.substr(0,pos), rhs = "";
    if (line.length() > (pos + _delim.length())) { rhs = line.substr(pos+_delim.length(),line.npos); }
    string rhsFull = rhs;
    pos = rhs.find(" ");
    rhs = rhs.substr(0,pos);
    vars.push_back(var); rhss.push_back(rhs); rhsFulls.push_back(rhsFull);
  }
}


bool BulkStateField::parseCommon(const string& var, const string& rhs, const string& rhsFull)
{
  if (var == _prefix + "type") { _type = rhs; }
  else if (var == _prefix + _symbol + "Vals") { _initVals.clear(); loadVectorFromInputFile(rhsFull,_initVals); }
  else if (var == _prefix + _symbol + "Depths") { _initDepths.clear(); loadVectorFromInputFile(rhsFull,_initDepths); }
  else if (var == _prefix + "eTest") { _eTest = atof(rhs.c_str()); }
  else if (var == _prefix + "eTestAmp") { _eTestAmp = atof(rhs.c_str()); }
  else if (var == _prefix + "eTestPeriod") { _eTestPeriod = atof(rhs.c_str()); }
  else { return false; }
  return true;
}


PetscErrorCode BulkStateField::setup()
{
  PetscErrorCode ierr = 0;
  if (_type != "transient" && _type != "constant") {
    PetscPrintf(PETSC_COMM_WORLD,"Error: %stype must be off, transient or constant (not %s).\n",_prefix.c_str(),_type.c_str());
    assert(0);
  }
  if (_initVals.size() != _initDepths.size()) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: %s%sVals and %s%sDepths must have the same length.\n",_prefix.c_str(),_symbol.c_str(),_prefix.c_str(),_symbol.c_str());
    assert(0);
  }
  if (_eTest >= 0 && !(_eTestPeriod > 0)) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: %seTestPeriod must be positive.\n",_prefix.c_str());
    assert(0);
  }
  ierr = VecDuplicate(_D->_y,&_state); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) _state, _symbol.c_str()); CHKERRQ(ierr);
  ierr = VecDuplicate(_D->_y,&_rate); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) _rate, (_symbol + "_t").c_str()); CHKERRQ(ierr);
  ierr = VecSet(_rate,0.0); CHKERRQ(ierr);
  ierr = VecDuplicate(_D->_y,&_e); CHKERRQ(ierr);
  ierr = VecSet(_e,0.0); CHKERRQ(ierr);

  // initial state: depth profile (or the default), then a file, then (on a restart) the checkpoint
  if (_initVals.size() > 0) { ierr = setVec(_state,_D->_z,_initVals,_initDepths); CHKERRQ(ierr); }
  else { ierr = initialDefault(_state); CHKERRQ(ierr); }
  ierr = loadVecFromInputFile(_state,_inputDir,_prefix + _symbol); CHKERRQ(ierr);
  if (_D->_restartFromChkpt) { ierr = loadCheckpoint(); CHKERRQ(ierr); }
  ierr = clamp(); CHKERRQ(ierr);
  return ierr;
}


PetscErrorCode BulkStateField::computeMaxTimeStep(const BulkInputs& in, PetscScalar& maxDeltaT)
{
  maxDeltaT = PETSC_MAX_REAL;
  return 0;
}


PetscErrorCode BulkStateField::writeContextExtra(PetscViewer& ascii, PetscViewer& viewer) { return 0; }
PetscErrorCode BulkStateField::clamp() { return 0; }
PetscErrorCode BulkStateField::initialDefault(Vec& state) { return VecSet(state,0.0); }


PetscErrorCode BulkStateField::drivingStrainRate(const BulkInputs& in)
{
  PetscErrorCode ierr = 0;
  if (_eTest >= 0) {
    const PetscScalar e = _eTest*(1.0 + _eTestAmp*sin(2.0*M_PI*in.time/_eTestPeriod));
    ierr = VecSet(_e,e); CHKERRQ(ierr);
  }
  else {
    const Vec rate = drivingRate(in);
    if (rate == NULL) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_NULL,"%s: its driving strain rate is not computed",_name.c_str()); }
    ierr = VecCopy(rate,_e); CHKERRQ(ierr);
    ierr = VecScale(_e,1e-3); CHKERRQ(ierr); // 1e-3/s -> 1/s
  }
  return ierr;
}


PetscErrorCode BulkStateField::initiateIntegrand(const PetscScalar time, map<string,Vec>& varEx)
{
  PetscErrorCode ierr = 0;
  if (_type != "transient") { return ierr; }
  if (varEx.find(_name) != varEx.end()) { ierr = VecCopy(_state,varEx[_name]); CHKERRQ(ierr); }
  else {
    Vec v;
    ierr = VecDuplicate(_state,&v); CHKERRQ(ierr);
    ierr = VecCopy(_state,v); CHKERRQ(ierr);
    varEx[_name] = v;
  }
  return ierr;
}


PetscErrorCode BulkStateField::updateFields(const PetscScalar time, const map<string,Vec>& varEx)
{
  PetscErrorCode ierr = 0;
  map<string,Vec>::const_iterator it = varEx.find(_name);
  if (it == varEx.end()) { return ierr; }
  ierr = VecCopy(it->second,_state); CHKERRQ(ierr);
  ierr = clamp(); CHKERRQ(ierr);
  return ierr;
}


PetscErrorCode BulkStateField::d_dt(const BulkInputs& in, map<string,Vec>& dvarEx)
{
  PetscErrorCode ierr = 0;
  if (_type != "transient") { return ierr; }
  ierr = drivingStrainRate(in); CHKERRQ(ierr);
  ierr = computeRate(in); CHKERRQ(ierr);
  ierr = VecCopy(_rate,dvarEx[_name]); CHKERRQ(ierr);
  return ierr;
}


// With a non-empty timeIntInds the state joins the step-size control with scale 1 (states are O(1)).
PetscErrorCode BulkStateField::addErrorControl(vector<string>& errInds, vector<double>& scale) const
{
  PetscErrorCode ierr = 0;
  if (_type != "transient" || errInds.empty()) { return ierr; }
  if (std::find(errInds.begin(),errInds.end(),_name) != errInds.end()) { return ierr; }
  while (scale.size() < errInds.size()) { scale.push_back(1.0); }
  errInds.push_back(_name);
  scale.push_back(1.0);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Note: %s added to timeIntInds (error scale 1).\n",_name.c_str()); CHKERRQ(ierr);
  return ierr;
}


PetscErrorCode BulkStateField::writeContext(PetscViewer& viewer)
{
  PetscErrorCode ierr = 0;
  PetscViewer ascii;
  string str = _outputDir + _name + ".txt";
  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&ascii); CHKERRQ(ierr);
  ierr = PetscViewerSetType(ascii,PETSCVIEWERASCII); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(ascii,FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(ascii,str.c_str()); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"%stype = %s\n",_prefix.c_str(),_type.c_str()); CHKERRQ(ierr);
  if (_eTest >= 0) {
    ierr = PetscViewerASCIIPrintf(ascii,"%seTest = %.15e # (1/s) test: prescribed driving strain rate\n",_prefix.c_str(),_eTest); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(ascii,"%seTestAmp = %.15e\n",_prefix.c_str(),_eTestAmp); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(ascii,"%seTestPeriod = %.15e # (s)\n",_prefix.c_str(),_eTestPeriod); CHKERRQ(ierr);
  }
  ierr = writeContextExtra(ascii,viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&ascii); CHKERRQ(ierr);
  return ierr;
}


PetscErrorCode BulkStateField::writeStep(PetscViewer& viewer)
{
  PetscErrorCode ierr = 0;
  ierr = PetscViewerHDF5PushGroup(viewer,("/" + _name).c_str()); CHKERRQ(ierr);
  ierr = PetscViewerHDF5PushTimestepping(viewer); CHKERRQ(ierr);
  ierr = VecView(_state,viewer); CHKERRQ(ierr);
  ierr = VecView(_rate,viewer); CHKERRQ(ierr);
  ierr = PetscViewerHDF5PopTimestepping(viewer); CHKERRQ(ierr);
  ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);
  return ierr;
}


PetscErrorCode BulkStateField::writeCheckpoint(PetscViewer& viewer)
{
  PetscErrorCode ierr = 0;
  ierr = PetscViewerHDF5PushGroup(viewer,("/" + _name).c_str()); CHKERRQ(ierr);
  ierr = VecView(_state,viewer); CHKERRQ(ierr);
  ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);
  return ierr;
}


PetscErrorCode BulkStateField::loadCheckpoint()
{
  PetscErrorCode ierr = 0;
  string fileName = _outputDir + "checkpoint.h5";
  PetscViewer viewer;
  ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,fileName.c_str(),FILE_MODE_READ,&viewer); CHKERRQ(ierr);
  ierr = PetscViewerHDF5PushGroup(viewer,("/" + _name).c_str()); CHKERRQ(ierr);
  ierr = VecLoad(_state,viewer); CHKERRQ(ierr);
  ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  return ierr;
}


PetscErrorCode refuseBulkStates(const char* file, const string& delim, const string& mediator)
{
  PetscErrorCode ierr = 0;
  const char* prefixes[] = {"hard_", "water_", "fabric_", "cement_", "seg_"};
  ifstream infile(file);
  string line;
  while (getline(infile, line)) {
    size_t pos = line.find(delim);
    if (pos == string::npos) { continue; }
    string var = line.substr(0,pos), rhs = line.substr(pos+delim.length());
    rhs = rhs.substr(0,rhs.find(" "));
    for (size_t k = 0; k < sizeof(prefixes)/sizeof(prefixes[0]); k++) {
      if (var == string(prefixes[k]) + "type" && rhs != "off") {
        PetscPrintf(PETSC_COMM_WORLD,"Error: %s = %s needs the power-law quasi-dynamic mediator; %s does not evolve it.\n",
          var.c_str(),rhs.c_str(),mediator.c_str());
        assert(0);
      }
    }
  }
  return ierr;
}
