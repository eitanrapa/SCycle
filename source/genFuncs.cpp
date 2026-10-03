#include "genFuncs.hpp"
#include <cstdio>
#include <cstring>
#include <signal.h>

using namespace std;

// check if file exists
bool doesFileExist(const string fileName)
{
    ifstream infile(fileName.c_str());
    return infile.good();
};


// clean up a C++ std library vector of PETSc Vecs
void destroyVector(vector<Vec>& vec)
{
  for(vector<Vec>::size_type i = 0; i != vec.capacity(); i++) {
    VecDestroy(&vec[i]);
  }
}

void destroyVector(map<string,Vec>& vec)
{
  for (map<string,Vec>::iterator it = vec.begin(); it!=vec.end(); it++ ) {
    VecDestroy(&vec[it->first]);
  }
}


// Print out a vector with 15 significant figures.
void printVec(Vec vec)
{
  //~ const char *name;
  //~ PetscObjectGetName((PetscObject) vec, &name);
  //~ PetscPrintf(PETSC_COMM_WORLD,"Vec %s: ",name);

  PetscInt Ii,Istart,Iend;
  PetscScalar v;
  VecGetOwnershipRange(vec,&Istart,&Iend);
  for (Ii = Istart;Ii < Iend;Ii++)
  {
    VecGetValues(vec,1,&Ii,&v);
    PetscPrintf(PETSC_COMM_WORLD,"%.15e\n",v);
  }
  PetscPrintf(PETSC_COMM_WORLD,"\n");
}

// Print out (vec1 - vec2) with 15 significant figures.
void printVecsDiff(Vec vec1,Vec vec2)
{
  PetscInt Ii,Istart,Iend;
  PetscScalar v1,v2,v;
  VecGetOwnershipRange(vec1,&Istart,&Iend);
  for (Ii = Istart;Ii < Iend;Ii++)
  {
    VecGetValues(vec1,1,&Ii,&v1);
    VecGetValues(vec2,1,&Ii,&v2);
    v = v1 - v2;
    PetscPrintf(PETSC_COMM_WORLD,"%.15e\n",v);
  }
}

// Print out (vec1 + vec2) with 15 significant figures.
void printVecsSum(Vec vec1,Vec vec2)
{
  PetscInt Ii,Istart,Iend;
  PetscScalar v1,v2,v;
  VecGetOwnershipRange(vec1,&Istart,&Iend);
  for (Ii = Istart;Ii < Iend;Ii++)
  {
    VecGetValues(vec1,1,&Ii,&v1);
    VecGetValues(vec2,1,&Ii,&v2);
    v = v1 + v2;
    PetscPrintf(PETSC_COMM_WORLD,"%.15e\n",v);
  }
}

// check that a Vec contains no inf or NaN values
double anyIsnan(const Vec& vec, string str)
{
  PetscErrorCode ierr = 0;

  const PetscScalar *val;
  PetscInt Istart, Iend;
  VecGetOwnershipRange(vec,&Istart,&Iend);
  VecGetArrayRead(vec,&val);
  PetscInt Jj = 0;
  for (PetscInt Ii=Istart;Ii<Iend;Ii++) {
    if (PetscIsNanReal(val[Jj]) || PetscIsInfReal(val[Jj])) {
      PetscPrintf(PETSC_COMM_WORLD,"%s",str.c_str());
      assert(!PetscIsNanReal(val[Jj]));
      assert(!PetscIsInfReal(val[Jj]));
    }
    Jj++;
  }
  VecRestoreArrayRead(vec,&val);

  return ierr;
}

// check that a Vec contains no inf or NaN values
double changeAnyIsnan(Vec& vec, string str, float newVal)
{
  PetscErrorCode ierr = 0;

  PetscScalar *val;
  PetscInt Istart, Iend;
  VecGetOwnershipRange(vec,&Istart,&Iend);
  VecGetArray(vec,&val);
  PetscInt Jj = 0;
  PetscInt hasBeenPrinted = 0;
  for (PetscInt Ii=Istart;Ii<Iend;Ii++) {
    if (PetscIsNanReal(val[Jj]) || PetscIsInfReal(val[Jj])) {
      if (hasBeenPrinted == 0) {
        PetscPrintf(PETSC_COMM_WORLD,"%s",str.c_str());
        hasBeenPrinted = 1;
      }
      val[Jj] = newVal;
      //~ assert(!PetscIsNanReal(val[Jj]));
      //~ assert(!PetscIsInfReal(val[Jj]));
    }
    Jj++;
  }
  VecRestoreArray(vec,&val);

  return ierr;
}

// write a single vector to file in HDF5 format
PetscErrorCode writeVec_hdf5(Vec vec, const string outFileName, const string group, const string objectName)
{
  PetscErrorCode err = 0;
  PetscViewer viewer;
  err = PetscViewerHDF5Open(PETSC_COMM_WORLD, outFileName.c_str(), FILE_MODE_APPEND, &viewer);CHKERRQ(err);
  err = PetscViewerHDF5PushGroup(viewer, group.c_str());CHKERRQ(err);
  err = PetscViewerHDF5PushTimestepping(viewer);CHKERRQ(err);
  PetscObjectSetName((PetscObject) vec, objectName.c_str());
  err = VecView(vec, viewer);CHKERRQ(err);
  err = PetscViewerHDF5PopTimestepping(viewer);CHKERRQ(err);
  err = PetscViewerHDF5PopGroup(viewer);CHKERRQ(err);
  PetscViewerDestroy(&viewer);

  return err;
}

// write a single vector to file in binary format
PetscErrorCode writeVec(Vec vec, const string filename)
{
  PetscErrorCode ierr = 0;
  PetscViewer    viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename.c_str(),FILE_MODE_WRITE,&viewer);
  ierr = VecView(vec,viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  return ierr;
}

// append a single vector to file in binary format
PetscErrorCode writeVecAppend(Vec vec,const string filename)
{
  PetscErrorCode ierr = 0;
  PetscViewer    viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename.c_str(),FILE_MODE_APPEND,&viewer);
  ierr = VecView(vec,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  return ierr;
}

// Write matrix to file in binary format
PetscErrorCode writeMat(Mat mat, string filename)
{
  PetscErrorCode ierr = 0;
  PetscViewer viewer;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename.c_str(),FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = MatView(mat,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  return ierr;
}

// initiate PetscViewer
PetscViewer initiateViewer(string filename)
{
  PetscViewer vw;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename.c_str(),FILE_MODE_WRITE,&vw);
  return vw;
}

// sappend PetscViewer(s)
PetscErrorCode appendViewer(PetscViewer& vw,const string filename)
{
  PetscErrorCode ierr = 0;
  ierr = PetscViewerDestroy(&vw); CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename.c_str(),FILE_MODE_APPEND,&vw);
  CHKERRQ(ierr);
  return ierr;
}

PetscErrorCode appendViewers(map<string,PetscViewer>& vwL,const string dir)
{
  PetscErrorCode ierr = 0;
  for (map<string,PetscViewer>::iterator it=vwL.begin(); it!=vwL.end(); it++ ) {
    ierr = PetscViewerDestroy(&vwL[it->first]); CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,(dir + it->first).c_str(),FILE_MODE_APPEND,&vwL[it->first]); CHKERRQ(ierr);
  }
  return ierr;
}


/* Print all entries of 2D DMDA global vector to stdout, including which
 * processor each entry lives on, and the corresponding subscripting indices
 */
PetscErrorCode printf_DM_2d(const Vec gvec, const DM dm)
{
  PetscErrorCode ierr = 0;
  #if VERBOSE > 2
    PetscPrintf(PETSC_COMM_WORLD,"Starting main::printf_DM_2d in fault.cpp.\n");
  #endif

  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscInt i,j,mStart,m,nStart,n; // for for loops below
  DMDAGetCorners(dm,&mStart,&nStart,0,&m,&n,0);

  PetscScalar **gxArr;
  DMDAVecGetArray(dm,gvec,&gxArr);
  for (j=nStart;j<nStart+n;j++) {
    for (i=mStart;i<mStart+m;i++) {
      PetscPrintf(PETSC_COMM_SELF,"%i: gxArr[%i][%i] = %g\n",
        rank,j,i,gxArr[j][i]);
    }
  }
  DMDAVecRestoreArray(dm,gvec,&gxArr);

  #if VERBOSE > 2
    PetscPrintf(PETSC_COMM_WORLD,"Ending main::printf_DM_2d in fault.cpp.\n");
  #endif
  return ierr;
}

// computes || vec1 - vec2||_mat (matrix norm of vec1 - vec2)
double computeNormDiff_Mat(const Mat& mat,const Vec& vec1,const Vec& vec2)
{
  PetscErrorCode ierr = 0;
  Vec diff;
  ierr = VecDuplicate(vec1,&diff);CHKERRQ(ierr);
  ierr = VecWAXPY(diff,-1.0,vec1,vec2);CHKERRQ(ierr);

  PetscScalar diffErr = computeNorm_Mat(mat,diff);
  PetscScalar vecErr = computeNorm_Mat(mat,vec1);
  PetscScalar err = diffErr/vecErr;

  VecDestroy(&diff);

  return err;
}

// computes || vec ||_mat
double computeNorm_Mat(const Mat& mat,const Vec& vec)
{
  PetscErrorCode ierr = 0;
  Vec Matxvec;
  ierr = VecDuplicate(vec,&Matxvec);CHKERRQ(ierr);
  ierr = MatMult(mat,vec,Matxvec);CHKERRQ(ierr);

  PetscScalar err;
  ierr = VecDot(vec,Matxvec,&err);CHKERRQ(ierr);

  VecDestroy(&Matxvec);
  return err;
}

// computes || vec1 - vec2 ||_2 / sqrt(length(vec1))
double computeNormDiff_2(const Vec& vec1,const Vec& vec2)
{
  PetscErrorCode ierr = 0;
  Vec diff;
  ierr = VecDuplicate(vec1,&diff);CHKERRQ(ierr);
  ierr = VecWAXPY(diff,-1.0,vec1,vec2);CHKERRQ(ierr);

  PetscScalar err;
  ierr = VecNorm(diff,NORM_2,&err);CHKERRQ(ierr);

  PetscInt len;
  ierr = VecGetSize(vec1,&len);CHKERRQ(ierr);
  err = err/sqrt(len);

  VecDestroy(&diff);

  return err;
}

// computes || vec1 - vec2 ||_2 / || vec1 ||_2
double computeNormDiff_L2_scaleL2(const Vec& vec1,const Vec& vec2)
{
  PetscErrorCode ierr = 0;
  Vec diff;
  ierr = VecDuplicate(vec1,&diff);CHKERRQ(ierr);
  ierr = VecWAXPY(diff,-1.0,vec1,vec2);CHKERRQ(ierr);

  PetscScalar err;
  ierr = VecNorm(diff,NORM_2,&err);CHKERRQ(ierr);

  PetscScalar len;
  ierr = VecNorm(vec1,NORM_2,&len);CHKERRQ(ierr);
  err = err/len;

  VecDestroy(&diff);

  return err;
}

// computes max( abs( vec1 - vec2 ./ / vec1 ))
double computeMaxDiff_scaleVec1(const Vec& vec1,const Vec& vec2)
{
  PetscErrorCode ierr = 0;
  PetscScalar err = 0;
  Vec diff;
  ierr = VecDuplicate(vec1,&diff);                                      CHKERRQ(ierr);
  ierr = VecWAXPY(diff,-1.0,vec1,vec2);                                 CHKERRQ(ierr);
  ierr = VecAbs(diff);                                                  CHKERRQ(ierr);
  ierr = VecPointwiseDivide(diff,diff,vec1);                            CHKERRQ(ierr);
  ierr = VecMax(diff,NULL,&err);                                        CHKERRQ(ierr);
  assert(!PetscIsInfReal(err));

  VecDestroy(&diff);

  return err;
}

// out = vecL' x A x vecR
double multVecMatsVec(const Vec& vecL, const Mat& A, const Vec& vecR)
{
  PetscErrorCode ierr = 0;
  double out = 0;
  Vec temp;
  VecDuplicate(vecL,&temp);
  ierr = MatMult(A,vecR,temp); CHKERRQ(ierr);
  ierr = VecDot(vecL,temp,&out); CHKERRQ(ierr);
  VecDestroy(&temp);

  return out;
}

// out = vecL' x A x B x vecR
double multVecMatsVec(const Vec& vecL, const Mat& A, const Mat& B, const Vec& vecR)
{
  PetscErrorCode ierr = 0;
  double out = 0;
  Mat AB;
  MatMatMult(A,B,MAT_INITIAL_MATRIX,1.0,&AB);

  Vec temp;
  VecDuplicate(vecL,&temp);

  ierr = MatMult(AB,vecR,temp); CHKERRQ(ierr);
  ierr = VecDot(vecL,temp,&out); CHKERRQ(ierr);

  VecDestroy(&temp);
  MatDestroy(&AB);
  return out;
}

// out = vecL' x A x B x C x vecR
double multVecMatsVec(const Vec& vecL, const Mat& A, const Mat& B, const Mat& C, const Vec& vecR)
{
  PetscErrorCode ierr = 0;
  double out = 0;

  Mat ABC;
  MatMatMatMult(A,B,C,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&ABC);

  Vec temp;
  VecDuplicate(vecL,&temp);
  ierr = MatMult(ABC,vecR,temp); CHKERRQ(ierr);
  ierr = VecDot(vecL,temp,&out); CHKERRQ(ierr);

  VecDestroy(&temp);
  MatDestroy(&ABC);
  return out;
}

// out = vecL' x A x B x C x vecR
double multVecMatsVec(const Vec& vecL, const Mat& A, const Mat& B, const Mat& C, const Mat& D, const Vec& vecR)
{
  PetscErrorCode ierr = 0;
  double out = 0;
  Mat ABC,ABCD;
  MatMatMatMult(A,B,C,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&ABC);
  MatMatMult(ABC,D,MAT_INITIAL_MATRIX,1.0,&ABCD);
  MatDestroy(&ABC);

  Vec temp;
  VecDuplicate(vecL,&temp);

  ierr = MatMult(ABCD,vecR,temp); CHKERRQ(ierr);
  ierr = VecDot(vecL,temp,&out); CHKERRQ(ierr);

  VecDestroy(&temp);
  MatDestroy(&ABCD);
  return out;
}

// out = A x B x vecR; assumes vecR is DIFFERENT from out
PetscErrorCode multMatsVec(Vec& out, const Mat& A, const Mat& B, const Vec& vecR)
{
  PetscErrorCode ierr = 0;
  Vec BvecR;
  VecDuplicate(vecR,&BvecR);
  ierr = MatMult(B,vecR,BvecR); CHKERRQ(ierr);
  ierr = MatMult(A,BvecR,out); CHKERRQ(ierr);
  ierr = VecDestroy(&BvecR);
  return ierr;
}

// out = A x B x vecR; assumes vecR is IDENTICAL to out
PetscErrorCode multMatsVec(const Mat& A, const Mat& B, Vec& vecR)
{
  PetscErrorCode ierr = 0;
  Vec BvecR;
  VecDuplicate(vecR,&BvecR);
  ierr = MatMult(B,vecR,BvecR); CHKERRQ(ierr);
  ierr = MatMult(A,BvecR,vecR); CHKERRQ(ierr);
  ierr = VecDestroy(&BvecR);
  return ierr;
}

// log10(out) = a*log10(vec1) + b*log10(vec2)
// out may not be the same as vec1 or vec2
PetscErrorCode MyVecLog10AXPBY(Vec& out,const double a, const Vec& vec1, const double b, const Vec& vec2)
{

  if (out == NULL) {
    VecDuplicate(vec1,&out);
  }

  // compute effective viscosity
  PetscScalar *outA;
  PetscScalar const *vec1A,*vec2A;
  PetscInt Ii,Istart,Iend,Jj = 0;
  VecGetOwnershipRange(vec1,&Istart,&Iend);
  VecGetArrayRead(vec1,&vec1A);
  VecGetArrayRead(vec2,&vec2A);
  VecGetArray(out,&outA);

  for (Ii = Istart;Ii < Iend;Ii++) {
   PetscScalar log10Out = a*log10(vec1A[Jj]) + b*log10(vec2A[Jj]);
    outA[Jj] = pow(10.,log10Out);
    Jj++;
  }
  VecRestoreArrayRead(vec1,&vec1A);
  VecRestoreArrayRead(vec2,&vec2A);
  VecRestoreArray(out,&outA);

  return 0;
}

// loads a PETSc Vec from a binary file
// Note: memory for out MUST be allocated before calling this function
PetscErrorCode loadVecFromInputFile(Vec& out,const string inputDir, const string fieldName)
{
  PetscErrorCode ierr = 0;
  bool fileExists = 0;
  ierr = loadVecFromInputFile(out,inputDir,fieldName,fileExists); CHKERRQ(ierr);
  return ierr;
}


// loads a PETSc Vec from a binary file
// Note: memory for out MUST be allocated before calling this function
PetscErrorCode loadVecFromInputFile(Vec& out,const string inputDir, const string fieldName, bool& fileExists)
{
  PetscErrorCode ierr = 0;
  #if VERBOSE > 1
    string funcName = "loadFieldsFromFiles";
    string fileName = "genFuncs.cpp";
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Starting %s in %s.\n",funcName.c_str(),fileName.c_str());CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  Attempting to load: %s%s\n",inputDir.c_str(),fieldName.c_str());CHKERRQ(ierr);
  #endif

  string vecSourceFile = inputDir + fieldName;

  fileExists = doesFileExist(vecSourceFile);
  if (fileExists) {
    PetscPrintf(PETSC_COMM_WORLD,"Note: Loading Vec from file: %s\n",vecSourceFile.c_str());
    PetscViewer inv;
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,vecSourceFile.c_str(),FILE_MODE_READ,&inv);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(inv,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
    ierr = VecLoad(out,inv);CHKERRQ(ierr);
    PetscViewerPopFormat(inv);
    PetscViewerDestroy(&inv);
  }
  else {
    PetscPrintf(PETSC_COMM_WORLD,"Warning: File not found: %s\n",vecSourceFile.c_str());
  }

  #if VERBOSE > 1
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Ending %s in %s.\n",funcName.c_str(),fileName.c_str());CHKERRQ(ierr);
  #endif
  return ierr;
}

// Split the contents of an input-file list "[a b c]" into whitespace-separated tokens.
// Anything after ']' is ignored; if there is no ']', parsing stops at a '#' comment.
static void splitInputFileList(const string& str,vector<string>& tokens)
{
  size_t first = str.find("[");
  size_t start = (first == string::npos) ? 0 : first + 1;
  size_t last = str.find("]",start);
  if (last == string::npos) { last = str.find("#",start); }
  string inner = str.substr(start, (last == string::npos) ? string::npos : last - start);
  istringstream iss(inner);
  string token;
  while (iss >> token) { tokens.push_back(token); }
}

// loads a std library vector from a list in the input file
PetscErrorCode loadVectorFromInputFile(const string& str,vector<double>& vec)
{
  PetscErrorCode ierr = 0;
  #if VERBOSE > 1
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Starting loadVectorFromInputFile in genFuncs.cpp.\n");CHKERRQ(ierr);
  #endif

  vector<string> tokens;
  splitInputFileList(str,tokens);
  for (size_t i = 0; i < tokens.size(); i++) { vec.push_back(atof(tokens[i].c_str())); }

  #if VERBOSE > 1
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Ending loadVectorFromInputFile in genFuncs.cpp.\n");CHKERRQ(ierr);
  #endif
  return ierr;
}

// load vector from input file
PetscErrorCode loadVectorFromInputFile(const string& str,vector<int>& vec)
{
  PetscErrorCode ierr = 0;
  #if VERBOSE > 1
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Starting loadVectorFromInputFile in genFuncs.cpp.\n");CHKERRQ(ierr);
  #endif

  vector<string> tokens;
  splitInputFileList(str,tokens);
  for (size_t i = 0; i < tokens.size(); i++) { vec.push_back(atoi(tokens[i].c_str())); }

  #if VERBOSE > 1
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Ending loadVectorFromInputFile in genFuncs.cpp.\n");CHKERRQ(ierr);
  #endif
  return ierr;
}

// loads a vector of strings
PetscErrorCode loadVectorFromInputFile(const string& str,vector<string>& vec)
{
  PetscErrorCode ierr = 0;
  #if VERBOSE > 1
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Starting loadVectorFromInputFile in genFuncs.cpp.\n");CHKERRQ(ierr);
  #endif

  splitInputFileList(str,vec);

  #if VERBOSE > 1
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Ending loadVectorFromInputFile in genFuncs.cpp.\n");CHKERRQ(ierr);
  #endif
  return ierr;
}

// creates a string containing the contents of C++ std library vector
string vector2str(const vector<double> vec)
{
  ostringstream ss;
  for (vector<double>::const_iterator Ii=vec.begin(); Ii != vec.end(); Ii++) {
    ss << " " << *Ii;
  }
  string str = "[" + ss.str() + "]";
  return str;
}

// creates a string containing the contents of C++ std library vector
string vector2str(const vector<int> vec)
{
  ostringstream ss;
  for (vector<int>::const_iterator Ii=vec.begin(); Ii != vec.end(); Ii++) {
    ss << " " << *Ii;
  }
  string str = "[" + ss.str() + "]";
  return str;
}

// creates a string containing the contents of C++ std library vector
string vector2str(const vector<string> vec)
{
  ostringstream ss;
  for (vector<string>::const_iterator Ii=vec.begin(); Ii != vec.end(); Ii++) {
    ss << " " << *Ii;
  }
  string str = "[" + ss.str() + "]";
  return str;
}

// prints an array to a single line in std out
PetscErrorCode printArray(const PetscScalar * arr,const PetscScalar len)
{
  PetscErrorCode ierr = 0;
  string funcName = "genFuncs::printArray";
  string fileName = "genFuncs.cpp";
  #if VERBOSE > 1
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Starting %s in %s\n",funcName.c_str(),fileName.c_str());
    CHKERRQ(ierr);
  #endif

  cout << "[";
  for(int i=0; i<len; i++) {
    cout << arr[i] << ",";
  }
  cout << "]" << endl;

  #if VERBOSE > 1
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Ending %s in %s\n",funcName.c_str(),fileName.c_str());
    CHKERRQ(ierr);
  #endif
  return ierr;
}

double MMS_test(const double y,const double z) {
  return y*10.0 + z;
}

double MMS_test(const double z) {
  return z;
}

// Fills vec with the linear interpolation between the pairs of points (vals,depths) in the direction set by input coord.
PetscErrorCode setVec(Vec& vec, const Vec& coord, vector<double>& vals,vector<double>& depths)
{
  PetscErrorCode ierr = 0;
  PetscInt       Istart,Iend,N;
  PetscScalar    v,z,z0,z1,v0,v1;

  // check the profile: as many values as depths, at least one, depths non-decreasing
  bool validProfile = vals.size() > 0 && vals.size() == depths.size();
  for (size_t ind = 1; validProfile && ind < depths.size(); ind++) {
    if (depths[ind] < depths[ind-1]) { validProfile = false; }
  }
  if (!validProfile) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: invalid <name>Vals/<name>Depths input profile (%zu values, %zu depths; "
      "need equally many, at least one, and non-decreasing depths). Values: %s Depths: %s\n",
      vals.size(),depths.size(),vector2str(vals).c_str(),vector2str(depths).c_str());
    assert(validProfile);
  }

  VecSet(vec,vals[0]);
  VecGetSize(coord,&N);
  // no interpolation to be done
  if (N == 1) {
    return ierr;
  }

  // build structure from generalized input
  size_t vecLen = depths.size();
  ierr = VecGetOwnershipRange(vec,&Istart,&Iend);CHKERRQ(ierr);
  PetscScalar *vecA;
  const PetscScalar *coordA;
  PetscInt Ii, Jj = 0;
  VecGetArray(vec,&vecA);
  VecGetArrayRead(coord,&coordA);
  for (Ii = Istart;Ii < Iend;Ii++)
  {
    z = coordA[Jj];
    for (size_t ind = 0; ind < vecLen - 1; ind++) {
      z0 = depths[0+ind];
      z1 = depths[0+ind+1];
      v0 = vals[0+ind];
      v1 = vals[0+ind+1];
      if (z1 == z0) { continue; } // repeated depth: a step in the profile, nothing to interpolate
      if (z >= z0 && z <= z1) {
        v = (v1 - v0)/(z1 - z0) * (z - z0) + v0;
        vecA[Jj] = v;
      }
      else if (z > z1 && ind == vecLen - 2) {
        v = (v1 - v0)/(z1 - z0) * (z - z0) + v0;
        vecA[Jj] = v;
      }
    }
    Jj++;
  }
  VecRestoreArray(vec,&vecA);
  VecRestoreArrayRead(coord,&coordA);
  return ierr;
}



PetscErrorCode mapToVec(Vec& vec, double(*func)(double,double), const Vec& yV, const double t)
{
  PetscErrorCode ierr = 0;
  PetscScalar const *y;
  PetscScalar *v;
  PetscInt Ii,Istart,Iend;
  ierr = VecGetOwnershipRange(vec,&Istart,&Iend); CHKERRQ(ierr);
  ierr = VecGetArrayRead(yV,&y);
  ierr = VecGetArray(vec,&v);
  PetscInt Jj = 0;
  for (Ii = Istart; Ii < Iend; Ii++) {
    v[Jj] = func(y[Jj],t);
    Jj++;
  }
  VecRestoreArrayRead(yV,&y);
  VecRestoreArray(vec,&v);
  return ierr;
}

PetscErrorCode mapToVec(Vec& vec, double(*func)(double),const Vec& yV)
{
  PetscErrorCode ierr = 0;
  PetscScalar const *y;
  PetscScalar *v;
  PetscInt Ii,Istart,Iend;
  ierr = VecGetOwnershipRange(vec,&Istart,&Iend); CHKERRQ(ierr);
  ierr = VecGetArrayRead(yV,&y);
  ierr = VecGetArray(vec,&v);
  PetscInt Jj = 0;
  for (Ii = Istart; Ii < Iend; Ii++) {
    v[Jj] = func(y[Jj]);
    Jj++;
  }
  VecRestoreArrayRead(yV,&y);
  VecRestoreArray(vec,&v);
  return ierr;
}

PetscErrorCode mapToVec(Vec& vec, double(*func)(double,double,double),
  const Vec& yV,const Vec& zV, const double t)
{
  PetscErrorCode ierr = 0;
  PetscScalar const *y,*z;
  PetscScalar *v;
  PetscInt Ii,Istart,Iend;
  ierr = VecGetOwnershipRange(vec,&Istart,&Iend); CHKERRQ(ierr);
  ierr = VecGetArrayRead(yV,&y);
  ierr = VecGetArrayRead(zV,&z);
  ierr = VecGetArray(vec,&v);
  PetscInt Jj = 0;
  for (Ii = Istart; Ii < Iend; Ii++) {
    v[Jj] = func(y[Jj],z[Jj],t);
    Jj++;
  }
  VecRestoreArrayRead(yV,&y);
  VecRestoreArrayRead(zV,&z);
  VecRestoreArray(vec,&v);
  return ierr;
}

PetscErrorCode mapToVec(Vec& vec, double(*func)(double,double),
  const Vec& yV,const Vec& zV)
{
  PetscErrorCode ierr = 0;
  PetscScalar const *y,*z;
  PetscScalar *v;
  PetscInt Ii,Istart,Iend;
  ierr = VecGetOwnershipRange(vec,&Istart,&Iend); CHKERRQ(ierr);
  ierr = VecGetArrayRead(yV,&y);
  ierr = VecGetArrayRead(zV,&z);
  ierr = VecGetArray(vec,&v);
  PetscInt Jj = 0;
  for (Ii = Istart; Ii < Iend; Ii++) {
    v[Jj] = func(y[Jj],z[Jj]);
    Jj++;
  }
  VecRestoreArrayRead(yV,&y);
  VecRestoreArrayRead(zV,&z);
  VecRestoreArray(vec,&v);
  return ierr;
}

// Map a function that acts on scalars to a 2D DMDA Vec
PetscErrorCode mapToVec(Vec& vec, double(*func)(double,double), const int N, const double dy, const double dz,DM da)
{
  // assumes vec has already been created and it's size has been allocated
  PetscErrorCode ierr = 0;

  PetscInt zS,yS,zn,yn;
  DMDAGetCorners(da, &zS, &yS, 0, &zn, &yn, 0);
  PetscInt zE = zS + zn;
  PetscInt yE = yS + yn;

  PetscScalar** arr;
  ierr = DMDAVecGetArray(da, vec, &arr);CHKERRQ(ierr);

  PetscInt yI,zI;
  PetscScalar y,z;
    for (yI = yS; yI < yE; yI++) {
      for (zI = zS; zI < zE; zI++) {
        y = yI * dy;
        z = zI * dz;
        arr[yI][zI] = func(y,z);
      }
    }

  ierr = DMDAVecRestoreArray(da, vec, &arr);CHKERRQ(ierr);
  VecAssemblyBegin(vec);
  VecAssemblyEnd(vec);

  return ierr;
}

// Map a function that acts on scalars to a 1DD DMDA Vec
PetscErrorCode mapToVec(Vec& vec, double(*func)(double), const int N, const double dz,DM da)
{
  // assumes vec has already been created and it's size has been allocated
  PetscErrorCode ierr = 0;

  PetscInt zS,zn;
  DMDAGetCorners(da, &zS, 0, 0, &zn, 0, 0);
  PetscInt zE = zS + zn;

  PetscScalar* arr;
  ierr = DMDAVecGetArray(da, vec, &arr);CHKERRQ(ierr);

  PetscInt zI;
  PetscScalar z;
  for (zI = zS; zI < zE; zI++) {
    z = zI * dz;
    arr[zI] = func(z);
  }

  ierr = DMDAVecRestoreArray(da, vec, &arr);CHKERRQ(ierr);
  VecAssemblyBegin(vec);
  VecAssemblyEnd(vec);

  return ierr;
}

// Map a function that acts on scalars to a 2D DMDA Vec
PetscErrorCode mapToVec(Vec& vec, double(*func)(double,double,double), const int N, const double dy, const double dz,const double t,DM da)
{
  // assumes vec has already been created and it's size has been allocated
  PetscErrorCode ierr = 0;

  PetscInt zS,yS,zn,yn;
  DMDAGetCorners(da, &zS, &yS, 0, &zn, &yn, 0);
  PetscInt zE = zS + zn;
  PetscInt yE = yS + yn;

  PetscScalar** arr;
  ierr = DMDAVecGetArray(da, vec, &arr);CHKERRQ(ierr);

  PetscInt yI,zI;
  PetscScalar y,z;
    for (yI = yS; yI < yE; yI++) {
      for (zI = zS; zI < zE; zI++) {
        y = yI * dy;
        z = zI * dz;
        arr[yI][zI] = func(y,z,t);
      }
    }

  ierr = DMDAVecRestoreArray(da, vec, &arr);CHKERRQ(ierr);
  VecAssemblyBegin(vec);
  VecAssemblyEnd(vec);

  return ierr;
}

// Print out a vector with 15 significant figures.
void printVec(const Vec vec,const DM da)
{
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscInt zS,yS,zn,yn;
  DMDAGetCorners(da, &zS, &yS, 0, &zn, &yn, 0);
  PetscInt zE = zS + zn;
  PetscInt yE = yS + yn;

  PetscScalar** arr;
  DMDAVecGetArray(da, vec, &arr);

  PetscInt yI,zI;
    for (yI = yS; yI < yE; yI++) {
      for (zI = zS; zI < zE; zI++) {
        PetscPrintf(PETSC_COMM_SELF,"%i: f(%i,%i) = %.2f\n",rank,yI,zI,arr[yI][zI]);
      }
    }

  DMDAVecRestoreArray(da, vec, &arr);
  VecAssemblyBegin(vec);
  VecAssemblyEnd(vec);
}

// repmat for vecs (i.e. vec -> [vec vec])
// Note: out must already be allocated onto processors
// n = # of repeats
PetscErrorCode repVec(Vec& out, const Vec& in, const PetscInt n)
{
  PetscErrorCode ierr = 0;
  PetscInt N,Istart,Iend,Ii;
  PetscScalar v = 0.0;
  vector<PetscScalar> vals(n);
  vector<PetscInt>    inds(n);

  VecGetSize(in,&N);
  VecGetOwnershipRange(in,&Istart,&Iend);
  for (Ii=Istart; Ii<Iend; Ii++ ) {
    ierr = VecGetValues(in,1,&Ii,&v);CHKERRQ(ierr);
    for (int i = 0; i < n; i++) {
      vals[i] = v;
      inds[i] = Ii + i*N;
    }
    ierr = VecSetValues(out,n,inds.data(),vals.data(),INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(out);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(out);CHKERRQ(ierr);

  return ierr;
}

// undoes repmat for vecs (i.e. [vec vec] -> vec)
// Note: out must already be allocated onto processors
// gIstart,gIend = global indices of in to store in out
PetscErrorCode sepVec(Vec& out, const Vec& in, const PetscInt gIstart, const PetscInt gIend)
{
  PetscErrorCode ierr = 0;
  PetscScalar v = 0.0;
  PetscInt Istart, Iend;

  VecGetOwnershipRange(in,&Istart,&Iend);
  for (PetscInt Ii = Istart; Ii < Iend; Ii++ ) {
    if (Ii >= gIstart && Ii < gIend) {
      ierr = VecGetValues(in,1,&Ii,&v);CHKERRQ(ierr);
      PetscInt Jj = Ii - gIstart;
      ierr = VecSetValue(out,Jj,v,INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  ierr = VecAssemblyBegin(out);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(out);CHKERRQ(ierr);

  return ierr;
}

// maps vec to bigger vec (i.e. vec -> [0 0 0, vec, 0 0 0])
// Note: out must already be allocated onto processors
// Istart,Iend = global indices of out to put the values of in into
PetscErrorCode distributeVec(Vec& out, const Vec& in, const PetscInt gIstart, const PetscInt gIend)
{
  PetscErrorCode ierr = 0;
  PetscScalar v = 0.0;
  PetscInt Istart, Iend, Ii, Jj;
  VecGetOwnershipRange(in,&Istart,&Iend);
  for (Ii=Istart; Ii<Iend; Ii++ ) {
    ierr = VecGetValues(in,1,&Ii,&v);CHKERRQ(ierr);
    Jj = Ii + gIstart;
    ierr = VecSetValue(out,Jj,v,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(out);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(out);CHKERRQ(ierr);

  return ierr;
}

//============================ Checkpoint Functions ===============================

// loading value from ASCII checkpoint file, which should only have one value (for scalar values: time and error)
PetscErrorCode loadValueFromCheckpoint(const string outputDir, const string filename, PetscScalar &value) {
  PetscErrorCode ierr = 0;
  string         checkpointFile = outputDir + filename;
  bool           fileExists = doesFileExist(checkpointFile);

  if (fileExists) {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Loading %s\n", filename.c_str()); CHKERRQ(ierr);
    PetscViewer viewer;
    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, checkpointFile.c_str()); CHKERRQ(ierr);
    PetscMPIInt    rank, size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (!rank) {
      ierr = PetscViewerASCIIRead(viewer, &value, 1, NULL, PETSC_SCALAR); CHKERRQ(ierr);
      MPI_Bcast(&value, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(&value, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    }
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  }
  else {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Warning: %s not found, setting to default value.\n", checkpointFile.c_str()); CHKERRQ(ierr);
  }

  return ierr;
}


// loads value from ASCII file (for integers)
PetscErrorCode loadValueFromCheckpoint(const string outputDir, const string filename, PetscInt &value) {
  PetscErrorCode ierr = 0;
  string         checkpointFile = outputDir + filename;
  bool           fileExists = doesFileExist(checkpointFile);

  if (fileExists) {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Loading %s\n", filename.c_str()); CHKERRQ(ierr);
    PetscViewer viewer;
    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, checkpointFile.c_str()); CHKERRQ(ierr);
    PetscMPIInt    rank, size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (!rank) {
      ierr = PetscViewerASCIIRead(viewer, &value, 1, NULL, PETSC_INT); CHKERRQ(ierr);
      MPI_Bcast(&value, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(&value, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    }
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  }
  else {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Warning: %s not found, setting to default value.\n", checkpointFile.c_str()); CHKERRQ(ierr);
  }

  return ierr;
}


//~ // initiate viewer to write and append vectors
//~ PetscErrorCode io_initiateWriteAppend(map<string, pair<PetscViewer,string>> &vwL, const string key, const Vec& vec, const string filename)
//~ {
  //~ PetscErrorCode ierr = 0;

  //~ // initiate viewer
  //~ PetscViewer viewer;
  //~ ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
  //~ vwL[key].first = viewer;
  //~ vwL[key].second = filename;
  //~ ierr = VecView(vec, viewer); CHKERRQ(ierr);
  //~ ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  //~ // reset to append mode
  //~ ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_APPEND, &vwL[key].first); CHKERRQ(ierr);
  //~ return ierr;
//~ }


//~ // append PetscVecs to existing files (saving new outputs to original data file during future checkpoints)
//~ PetscErrorCode initiate_appendVecToOutput(map<string, pair<PetscViewer, string>> &vwL, const string key, const Vec &vec, const string filename) {
  //~ PetscErrorCode ierr = 0;

  //~ // initiate viewer
  //~ PetscViewerCreate(PETSC_COMM_WORLD, &vwL[key].first);
  //~ vwL[key].second = filename;
  //~ PetscViewerSetType(vwL[key].first, PETSCVIEWERBINARY);
  //~ PetscViewerFileSetMode(vwL[key].first, FILE_MODE_APPEND);
  //~ PetscViewerFileSetName(vwL[key].first, filename.c_str());
  //~ VecView(vec,vwL[key].first);
  //~ return ierr;
//~ }

PetscErrorCode writeASCII(const string outputDir, const string filename, PetscInt var,const string format) {
  PetscErrorCode ierr = 0;
  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, (outputDir + filename).c_str(), &viewer);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer, format.c_str(), var);CHKERRQ(ierr);
  PetscViewerDestroy(&viewer);
  return ierr;
}
PetscErrorCode writeASCII(const string outputDir, const string filename, PetscScalar var,const string format) {
  PetscErrorCode ierr = 0;
  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, (outputDir + filename).c_str(), &viewer);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer, format.c_str(), var);CHKERRQ(ierr);
  PetscViewerDestroy(&viewer);
  return ierr;
}

// initiate viewer and write to an ASCII file using the specified format
// mode can be FILE_MODE_WRITE or FILE_MODE_APPEND
// write ASCII file for scalar
PetscErrorCode initiateWriteASCII(const string outputDir, const string filename, const PetscFileMode mode, PetscViewer &viewer, const string format, PetscScalar var)
{
  PetscErrorCode ierr = 0;
  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);                   CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);                  CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, mode);                          CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, (outputDir + filename).c_str());CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer, format.c_str(), var);CHKERRQ(ierr);

  // ensure that viewer mode switches to append if it isn't that already
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_APPEND);              CHKERRQ(ierr);

  return ierr;
}


// write ASCII file for integer
PetscErrorCode initiateWriteASCII(const string outputDir, const string filename, const PetscFileMode mode, PetscViewer &viewer, const string format, PetscInt var)
{
  PetscErrorCode ierr = 0;
  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);                   CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);                  CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, mode);                          CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, (outputDir + filename).c_str());CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer, format.c_str(), var);CHKERRQ(ierr);

  // ensure that viewer mode switches to append if it isn't that already
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_APPEND);              CHKERRQ(ierr);

  return ierr;
}


// append PetscVecs to existing files (saving new outputs to original data file during future checkpoints)
PetscErrorCode initiate_writeVec_hdf5(map<string, pair<PetscViewer, string>> &vwL, const string key, const Vec &vec, const string filename, const PetscFileMode mode)
{
  PetscErrorCode ierr = 0;

  // initiate viewer
  ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(), mode, &vwL[key].first);CHKERRQ(ierr);
  ierr = PetscViewerHDF5PushTimestepping(vwL[key].first);             CHKERRQ(ierr);

  // write Vec for the first time
  ierr = VecView(vec,vwL[key].first);                                   CHKERRQ(ierr);
  ierr = PetscViewerHDF5PopTimestepping(vwL[key].first);              CHKERRQ(ierr);

  // ensure that viewer mode switches to append if it isn't that already
  ierr = PetscViewerDestroy(&vwL[key].first);
  ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_APPEND, &vwL[key].first);CHKERRQ(ierr);
  ierr = PetscViewerHDF5PushTimestepping(vwL[key].first);             CHKERRQ(ierr);

  return ierr;
}

// append PetscVecs to existing files (saving new outputs to original data file during future checkpoints)
PetscErrorCode initiate_appendVecToOutput(map<string, pair<PetscViewer, string>> &vwL, const string key, const Vec &vec, const string filename, const PetscFileMode mode)
{
  PetscErrorCode ierr = 0;

  // initiate viewer
  PetscViewerCreate(PETSC_COMM_WORLD, &vwL[key].first);
  vwL[key].second = filename;
  PetscViewerSetType(vwL[key].first, PETSCVIEWERBINARY);
  PetscViewerFileSetName(vwL[key].first, filename.c_str());
  PetscViewerFileSetMode(vwL[key].first, mode);

  // write Vec for the first time
  VecView(vec,vwL[key].first);

  // ensure that viewer mode switches to append if it isn't that already
  PetscViewerDestroy(&vwL[key].first);
  PetscViewerCreate(PETSC_COMM_WORLD, &vwL[key].first);
  PetscViewerSetType(vwL[key].first, PETSCVIEWERBINARY);
  PetscViewerFileSetName(vwL[key].first, filename.c_str());
  PetscViewerFileSetMode(vwL[key].first, FILE_MODE_APPEND);
  //~ PetscViewerFileSetMode(vwL[key].first, FILE_MODE_APPEND);

  return ierr;
}


// initiate viewer to write and append vectors
PetscErrorCode io_initiateWriteAppend(map<string, pair<PetscViewer,string>> &vwL, const string key, const Vec& vec, const string filename)
{
  PetscErrorCode ierr = 0;

  // initiate viewer
  PetscViewer viewer;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
  vwL[key].first = viewer;
  vwL[key].second = filename;
  ierr = VecView(vec, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  // reset to append mode
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_APPEND, &vwL[key].first); CHKERRQ(ierr);
  return ierr;
}






// max over entries of |a - b| / |b|, using |a - b| where b == 0
PetscReal maxRelativeDiff(const Vec& a, const Vec& b)
{
  const PetscScalar *aA,*bA;
  PetscInt Istart,Iend;
  VecGetOwnershipRange(b,&Istart,&Iend);
  VecGetArrayRead(a,&aA);
  VecGetArrayRead(b,&bA);
  PetscReal localMax = 0;
  for (PetscInt Jj = 0; Jj < Iend-Istart; Jj++) {
    PetscReal d = PetscAbsScalar(aA[Jj] - bA[Jj]);
    PetscReal scale = PetscAbsScalar(bA[Jj]);
    if (scale > 0) { d /= scale; }
    localMax = PetscMax(localMax,d);
  }
  VecRestoreArrayRead(a,&aA);
  VecRestoreArrayRead(b,&bA);
  PetscReal globalMax = 0;
  MPI_Allreduce(&localMax,&globalMax,1,MPIU_REAL,MPIU_MAX,PetscObjectComm((PetscObject)b));
  return globalMax;
}


// make an open HDF5 file consistent on disk (PetscViewerFlush does nothing for HDF5 viewers)
PetscErrorCode flushHDF5Viewer(PetscViewer viewer)
{
  PetscErrorCode ierr = 0;
  if (viewer == NULL) { return ierr; }
  hid_t fileId;
  ierr = PetscViewerHDF5GetFileId(viewer, &fileId); CHKERRQ(ierr);
  if (H5Fflush(fileId, H5F_SCOPE_GLOBAL) < 0) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_WRITE, "H5Fflush failed");
  }
  return ierr;
}

std::string checkpointTmpName(const std::string& outputDir)
{
  return outputDir + "checkpoint.h5.tmp";
}

// close the checkpoint written to the temporary file and move it over <outputDir>checkpoint.h5
PetscErrorCode commitCheckpoint(PetscViewer& viewer, const std::string& outputDir)
{
  PetscErrorCode ierr = 0;
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr); // collective: every rank has closed the file
  PetscMPIInt rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  int failed = 0;
  if (rank == 0) {
    const std::string tmp = checkpointTmpName(outputDir), final = outputDir + "checkpoint.h5";
    failed = (rename(tmp.c_str(), final.c_str()) != 0);
  }
  ierr = MPI_Bcast(&failed, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
  if (failed) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_WRITE, "could not rename checkpoint.h5.tmp to checkpoint.h5");
  }
  return ierr;
}


static volatile sig_atomic_t stopSignalReceived = 0; // set by the handler on this rank
static int stopSignalDecided = 0;                     // agreed by all ranks in checkStopRequest

extern "C" void stopSignalHandler(int sig)
{
  if (stopSignalReceived == sig) { signal(sig, SIG_DFL); raise(sig); } // asked twice: terminate now
  stopSignalReceived = sig;
}

void installStopHandler()
{
  struct sigaction action;
  memset(&action, 0, sizeof(action));
  action.sa_handler = stopSignalHandler;
  sigemptyset(&action.sa_mask);
  action.sa_flags = SA_RESTART; // system calls interrupted by the signal carry on
  sigaction(SIGTERM, &action, NULL);
  sigaction(SIGINT, &action, NULL);
}

bool checkStopRequest()
{
  if (stopSignalDecided != 0) { return true; }
  int local = stopSignalReceived, global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
  stopSignalDecided = global;
  return global != 0;
}

int stopSignal() { return stopSignalDecided; }

PetscErrorCode reportStop(PetscInt stepCount, PetscScalar time, bool checkpointed)
{
  PetscErrorCode ierr = 0;
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Stopping on signal %i at step %i, t = %.15e s: %s\n", stopSignalDecided,
                     stepCount, time, checkpointed ? "checkpoint written, run again to continue" : "no checkpoint written");
  CHKERRQ(ierr);
  return ierr;
}
