/*=========================================================================

  Program:   Visualization Toolkit
  Module:    VonMisesStress.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __VonMisesStress_h
#define __VonMisesStress_h

#include "vtkDataSetAlgorithm.h"

class VonMisesStress : public vtkDataSetAlgorithm
{
  
  public:
   
    vtkTypeMacro(VonMisesStress, vtkDataSetAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static VonMisesStress *New();
   

  protected:
    VonMisesStress();
    ~VonMisesStress();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

  private:
    VonMisesStress(const VonMisesStress&);  // Not implemented.
    void operator=(const VonMisesStress&);  // Not implemented.
  
};

#endif
