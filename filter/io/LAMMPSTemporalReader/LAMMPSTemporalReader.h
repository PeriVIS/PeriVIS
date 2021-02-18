/*
MIT License

Copyright (c) 2021 Michael Bussler

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#ifndef __LAMMPSTemporalReader_h
#define __LAMMPSTemporalReader_h

#include <vtkPolyDataAlgorithm.h>
#include <vtkSmartPointer.h>
#include <vector>
#include <string>

class vtkPolyData;
class vtkImplicitFunction;

#if defined(_MSC_VER) // visual studio is.. different!
 typedef          __int64 int64_t;
 typedef unsigned __int64 uint64_t;
#else
 #include <inttypes.h>
#endif

class Atom {
public:
	int id;
	int type;
	float pos[3];
	float nDataFields;
	float* data;
	Atom() : data(0) {};
	~Atom(){ if( data ) delete[] data; }
};

class Timestep {
public:
	int time;
	int nAtoms, nBonds;
	float box[6];
	std::vector<std::string> dataFieldNames;
	int nDataFields;
	int posField[3];
	int bondField[2];
	uint64_t atomStartLine;
	uint64_t bondsStartLine;
};

class LAMMPSTemporalReader : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(LAMMPSTemporalReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  static LAMMPSTemporalReader *New();
 
  // Description:
  // Specify file name of the .dat file.
  //vtkSetStringMacro(FileName);
  //vtkGetStringMacro(FileName);
  //vtkSetStringMacro(BondsFileName);
  //vtkGetStringMacro(BondsFileName);

  vtkGetMacro(CalculateBrokenBonds, bool);
  vtkSetMacro(CalculateBrokenBonds, bool);
  vtkGetMacro(CalculateBondLength, bool);
  vtkSetMacro(CalculateBondLength, bool);
  vtkGetMacro(CalculateRefBondLength, bool);
  vtkSetMacro(CalculateRefBondLength, bool);
  vtkGetMacro(LoadBonds, bool);
  vtkSetMacro(LoadBonds, bool);
  vtkGetMacro(MaxBondLength, double);
  vtkSetMacro(MaxBondLength, double);
  vtkGetMacro(FilterBondsByLength, bool);
  vtkSetMacro(FilterBondsByLength, bool);

  // Description:
  // Adds names of files to be read. The files are read in the order
  // they are added.
  virtual void AddFileName(const char* fname);

  // Description:
  // Remove all file names.
  virtual void RemoveAllFileNames();


protected:
  LAMMPSTemporalReader(vtkImplicitFunction *cf=NULL);
  ~LAMMPSTemporalReader();
 
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
  LAMMPSTemporalReader(const LAMMPSTemporalReader&);  // Not implemented.
  void operator=(const LAMMPSTemporalReader&);  // Not implemented.

  void ParseLAMMPSFile();
  void ParseBondsFile();
 
  //char* FileName;
  //char* BondsFileName;

  int NumberOfTimeSteps;
  std::vector<Timestep> m_timesteps;
  std::vector<double> TimestepValues;

  int FindClosestTimeStep(double requestedTimeValue);
  void ReadAtoms( int id, vtkPolyData* pd );
  void tokenize(const std::string line, std::vector<std::string> &tokens);
  void ReadBonds(int timestepToLoad, vtkPolyData * pd);

  vtkSmartPointer<vtkIdList> getConnectedPoints(vtkPolyData * pd, vtkIdType ptId, vtkSmartPointer<vtkIdList> bondIds = NULL);

  void WriteHeaderFile();
  bool ReadHeaderFile();
  int NumberOfComponents;
  bool bRelativePath;
  bool CalculateBrokenBonds;
  bool CalculateBondLength;
  bool CalculateRefBondLength;
  bool LoadBonds;
  bool FilterBondsByLength;
  double MaxBondLength;

  vtkSmartPointer<vtkPolyData> m_initialMesh;
  std::vector<std::string> m_FileNames;
  std::vector<std::string> m_BondFileNames;
};
 
#endif//__LAMMPSTemporalReader_h
