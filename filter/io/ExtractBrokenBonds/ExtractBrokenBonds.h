
#ifndef __ClipPolyData_h
#define __ClipPolyData_h

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkPolyData;
class vtkIdList;

class ExtractBrokenBonds : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(ExtractBrokenBonds,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  static ExtractBrokenBonds *New();

  vtkSetClampMacro( Samples,int,1,20);
  vtkGetMacro( Samples,int);

  vtkSetClampMacro( DiskSize, double, 0.0, 1.0);
  vtkGetMacro( DiskSize, double);

  vtkSetMacro( ScaleByBondlength, bool );
  vtkGetMacro( ScaleByBondlength, bool );

protected:
  ExtractBrokenBonds();
  ~ExtractBrokenBonds();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:

  vtkSmartPointer<vtkIdList> getConnectedPoints( vtkPolyData * pd, vtkIdType ptId, vtkSmartPointer<vtkIdList> bondIds = NULL);

  ExtractBrokenBonds(const ExtractBrokenBonds&);  // Not implemented.
  void operator=(const ExtractBrokenBonds&);  // Not implemented.

  int Samples;
  double DiskSize;
  bool ScaleByBondlength;
};

#endif
