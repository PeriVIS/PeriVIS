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


#ifndef __PeriConnectedComponents_h
#define __PeriConnectedComponents_h

#include <vtkPolyDataAlgorithm.h>
#include <vtkSmartPointer.h>

class vtkIdList;
class vtkPolyData;

class PeriConnectedComponents : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(PeriConnectedComponents,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  static PeriConnectedComponents *New();
 
  // Description:
  vtkSetClampMacro( MaxDamage,double,0.0,1.0);
  vtkGetMacro( MaxDamage,double);

  vtkSetMacro( MaxBondlength, double);
  vtkGetMacro( MaxBondlength, double);


protected:
  PeriConnectedComponents();
  ~PeriConnectedComponents();
 
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int FillInputPortInformation(int port, vtkInformation *info);

private:
  PeriConnectedComponents(const PeriConnectedComponents&);  // Not implemented.
  void operator=(const PeriConnectedComponents&);  // Not implemented.

  vtkSmartPointer<vtkIdList> 
    getConnectedPoints( vtkPolyData * pd, vtkIdType ptId, vtkSmartPointer<vtkIdList> bondIds = NULL);

  double MaxDamage;
  double MaxBondlength;
};
 
#endif//__PeriConnectedComponents_h
