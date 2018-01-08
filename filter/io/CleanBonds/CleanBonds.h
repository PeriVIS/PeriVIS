/*=========================================================================

Program:   Visualization Toolkit
Module:    CleanBonds.h

Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __CleanBonds_h
#define __CleanBonds_h

#include "vtkPolyDataAlgorithm.h"

class CleanBonds : public vtkPolyDataAlgorithm
{

public:

    vtkTypeMacro(CleanBonds, vtkPolyDataAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);

    // Description
    // Calculate Eigenvectors for tensor data
    static CleanBonds *New();


    vtkSetMacro( MaxBondLength, double);
    vtkGetMacro( MaxBondLength, double);

protected:
    CleanBonds();
    ~CleanBonds();

    //     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int FillInputPortInformation(int port, vtkInformation *info);

private:
    CleanBonds(const CleanBonds&);  // Not implemented.
    void operator=(const CleanBonds&);  // Not implemented.

    double MaxBondLength;
};

#endif
