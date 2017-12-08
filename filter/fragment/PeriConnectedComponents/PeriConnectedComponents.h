/*=========================================================================

Program:   Visualization Toolkit
Module:    PeriConnectedComponents.h

Copyright (c) Michael Buﬂler

=========================================================================*/

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
