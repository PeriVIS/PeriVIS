/*=========================================================================

Program:   Visualization Toolkit
Module:    CleanBonds.h

Copyright (c) Michael Bußler

=========================================================================*/

#include "CleanBonds.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkDataSet.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPolyData.h"
#include "vtkTransform.h"
#include "vtkSmartPointer.h"

vtkStandardNewMacro(CleanBonds);

//==============================================================================
CleanBonds::CleanBonds()
{
    this->MaxBondLength = 1.0;
}
//==============================================================================
CleanBonds::~CleanBonds()
{
    // clean up
}

//==============================================================================
int CleanBonds::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and output
    vtkDataSet *input = vtkDataSet::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output = vtkPolyData::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkPoints *newPoints;
    vtkCellArray *newLines, *newVerts;
    vtkIdType estimatedSize, numCells=input->GetNumberOfCells();
    vtkIdType numPts=input->GetNumberOfPoints();
    //vtkPoints *inPts=input->GetPoints();
    int numberOfPoints;
    vtkPointData *inPD=input->GetPointData(), *outPD = output->GetPointData();
    vtkCellData *inCD=input->GetCellData(), *outCD = output->GetCellData();

    if ( numPts < 1 /*|| inPts == NULL*/ )
    {
        vtkDebugMacro(<<"No data to convert");
        return 1;
    }

    // copy point data arrays
    int numArrays = inPD->GetNumberOfArrays();
    for( int i=0; i<numArrays; i++ ){
        int type = inPD->GetArray(i)->GetDataType();
        vtkDataArray *arr = vtkDataArray::CreateDataArray(type);
        arr->SetNumberOfComponents(inPD->GetArray(i)->GetNumberOfComponents());
        arr->SetName(inPD->GetArray(i)->GetName());
        outPD->AddArray(arr);
        arr->Delete();
    }

    // copy cell data arrays
    int numCellArrays = inCD->GetNumberOfArrays();
    for( int i=0; i<numCellArrays; i++ ){
        int type = inCD->GetArray(i)->GetDataType();
        vtkDataArray *arr = vtkDataArray::CreateDataArray(type);
        arr->SetNumberOfComponents(inCD->GetArray(i)->GetNumberOfComponents());
        arr->SetName(inCD->GetArray(i)->GetName());
        outCD->AddArray(arr);
        arr->Delete();
    }


    // copy all points and point data
    newPoints =  vtkPoints::New();

    vtkIdType vertexId = 0;

    for( vtkIdType ptId=0; ptId<numPts; ptId++) 
    {
        newPoints->InsertNextPoint( input->GetPoint(ptId));

        // copy point data
        for( int i=0; i<numArrays; i++ ){
            double data[9];
            inPD->GetArray(i)->GetTuple(ptId, data);
            outPD->GetArray(i)->InsertNextTuple(data);
        }

        this->UpdateProgress(0.5*ptId/numPts);
    }

    // currently only lines (bonds) are supported
    newLines = vtkCellArray::New();

    // mark connected points
    bool *connected = new bool [numPts];
    for( int i=0; i<numPts; i++)
        connected[i]=false;

    // maximum length of bonds to copy
    double maxLength2 = this->MaxBondLength*this->MaxBondLength;

    for( vtkIdType cellId=0; cellId<numCells; cellId++) 
    {
        vtkCell* cell = input->GetCell(cellId);
        if( cell->GetCellType() == VTK_LINE )
        {
            int pid1 = cell->GetPointId(0);
            int pid2 = cell->GetPointId(1);

            double p1[3], p2[3];
            input->GetPoint( pid1, p1);
            input->GetPoint( pid2, p2);

            double dist2 = vtkMath::Distance2BetweenPoints(p1, p2);

            if( dist2 < maxLength2 )
            {
                newLines->InsertNextCell(2);
                newLines->InsertCellPoint( cell->GetPointId(0) );
                newLines->InsertCellPoint( cell->GetPointId(1) );

                // copy cell data
                for( int i=0; i<numCellArrays; i++ ){
                    double data[9];
                    inCD->GetArray(i)->GetTuple(cellId, data);
                    outCD->GetArray(i)->InsertNextTuple(data);
                }

                connected[pid1] = true;
                connected[pid2] = true;
            }
        }
    }

    // generate vertices for unconnected atoms
    output->BuildLinks();

    //newVerts = vtkCellArray::New();

    //for( vtkIdType ptId=0; ptId<numPts; ptId++) 
    //{
    //    if (!connected[ptId]) // lone points
    //    {
    //        vtkIdType cellId = newVerts->InsertNextCell(1);
    //        newVerts->InsertCellPoint(ptId);

    //        // copy cell data
    //        for( int i=0; i<numCellArrays; i++ ){
    //            double data[9] = {0,0,0,0,0,0,0,0,0};
    //            outCD->GetArray(i)->InsertNextTuple(data);
    //        }
    //    }
    //}

    //
    // Update output and release memory
    //

    delete[] connected;

    output->SetPoints(newPoints);
    newPoints->Delete();

    output->SetLines(newLines);
    newLines->Delete();

    //output->SetVerts(newVerts);
    //newVerts->Delete();

    output->Squeeze();

    return 1;
}

//----------------------------------------------------------------------------
int CleanBonds::FillInputPortInformation(int, vtkInformation *info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
}

//==============================================================================
void CleanBonds::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}
