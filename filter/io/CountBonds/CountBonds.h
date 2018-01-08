
#ifndef __CountBonds_h
#define __CountBonds_h

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"

class CountBonds : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(CountBonds,vtkPolyDataAlgorithm);

  // Description:
  // Construct with user-specified implicit function; InsideOut turned off;
  // value set to 0.0; and generate clip scalars turned off.
  static CountBonds *New();

protected:
  CountBonds();
  ~CountBonds();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:

  vtkSmartPointer<vtkIdList> getConnectedPoints( vtkPolyData * pd, vtkIdType ptId, vtkSmartPointer<vtkIdList> bondIds = NULL);

  CountBonds(const CountBonds&);  // Not implemented.
  void operator=(const CountBonds&);  // Not implemented.
};

#endif
