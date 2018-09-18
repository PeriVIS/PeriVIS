/*=========================================================================

  Program:   Visualization Toolkit
  Module:    MakeTensor.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __MakeTensor_h
#define __MakeTensor_h

#include "vtkDataSetAlgorithm.h"

class MakeTensor : public vtkDataSetAlgorithm
{
  
  public:
   
    vtkTypeMacro(MakeTensor, vtkDataSetAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static MakeTensor *New();
   

  protected:
    MakeTensor();
    ~MakeTensor();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

  private:
    MakeTensor(const MakeTensor&);  // Not implemented.
    void operator=(const MakeTensor&);  // Not implemented.
  
};

#endif
