#include "ExtractBrokenBonds.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkGenericCell.h"
#include "vtkImplicitFunction.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkLine.h"
#include "vtkMergePoints.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkTriangle.h"
#include "vtkIncrementalPointLocator.h"
#include "vtkTransform.h"

#include <math.h>
#include "vtkSmartPointer.h"
#include <algorithm>

vtkStandardNewMacro(ExtractBrokenBonds);

//----------------------------------------------------------------------------
ExtractBrokenBonds::ExtractBrokenBonds()
{
    // by default, process active cell scalars
    this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS,
        vtkDataSetAttributes::SCALARS);

    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);

    this->Samples = 1;
    this->DiskSize = 0.1;
}

//----------------------------------------------------------------------------
ExtractBrokenBonds::~ExtractBrokenBonds()
{
}

//----------------------------------------------------------------------------
//
// Clip through data generating surface.
//
int ExtractBrokenBonds::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and output
    vtkPolyData *input = vtkPolyData::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output = vtkPolyData::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkIdType cellId, i, updateTime;
    vtkPoints *cellPts;
    vtkDataArray *clipScalars;
    vtkFloatArray *cellScalars;
    vtkCellArray *newTris, *newLines, *newPolys, *connList=NULL;
    vtkPoints *newPoints;
    double s;
    vtkIdType estimatedSize, numCells=input->GetNumberOfCells();
    vtkIdType numPts=input->GetNumberOfPoints();
    vtkPoints *inPts=input->GetPoints();
    int numberOfPoints;
    vtkPointData *inPD=input->GetPointData(), *outPD = output->GetPointData();
    vtkCellData *inCD=input->GetCellData(), *outCD = output->GetCellData();
    vtkCellData *outClippedCD = NULL;

    vtkDebugMacro(<< "Extracting broken bonds");

    vtkDataArray *bondActive = this->GetInputArrayToProcess(0, inputVector);

    // Initialize self; create output objects
    //
    if ( bondActive == NULL || numPts < 1 || inPts == NULL )
    {
        vtkDebugMacro(<<"No data.");
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

    // new point data
    newPoints =  vtkPoints::New();
    newTris = vtkCellArray::New();

    input->BuildLinks();

    // iterate over all bonds, insert a vertex for each broken bond
    for( cellId=0; cellId<numCells; cellId++) 
    {
        vtkCell* cell = input->GetCell(cellId);

        if( cell && cell->GetCellType() == VTK_LINE ) {

            double active = bondActive->GetTuple1(cellId);
            if( active < 1.0 ) { // 0.0 -> broken
                vtkSmartPointer<vtkIdList> pointIds = cell->GetPointIds();
                vtkIdType ptId1 = pointIds->GetId(0), 
                    ptId2 = pointIds->GetId(1);
                double pos1[3], pos2[3], dir[3], up[3], center[3], newPoint[3];
                vtkIdType centerPtId, prevPtId, currPtId, firstPtId;

                input->GetPoint( ptId1, pos1 );
                input->GetPoint( ptId2, pos2 );
                vtkMath::Subtract(pos2, pos1, dir);
                vtkMath::MultiplyScalar(dir, 0.5);
                vtkMath::Add(pos1, dir, center);

                centerPtId = newPoints->InsertNextPoint(center);

                up[0] = dir[1] - dir[2];
                up[1] = -dir[0];
                up[2] = dir[0];

                vtkMath::Normalize(up);
                
                if( this->ScaleByBondlength )
                    vtkMath::MultiplyScalar(up, vtkMath::Norm(dir)*2.0*this->DiskSize);
                else
                    vtkMath::MultiplyScalar(up, this->DiskSize);

                vtkMath::Add(center, up, newPoint);
                firstPtId = newPoints->InsertNextPoint(newPoint);
                prevPtId = firstPtId;

                vtkSmartPointer<vtkTransform> rot = vtkSmartPointer<vtkTransform>::New();
                rot->RotateWXYZ(360.0/this->Samples, dir);
                
                for( int i=0; i<this->Samples-1; i++)
                {
                    rot->TransformPoint(up, up);
                    vtkMath::Add(center, up, newPoint);
                    currPtId = newPoints->InsertNextPoint(newPoint);
                   
                    newTris->InsertNextCell(3);
                    newTris->InsertCellPoint(centerPtId);
                    newTris->InsertCellPoint(prevPtId);
                    newTris->InsertCellPoint(currPtId);

                    prevPtId = currPtId;
                }

                newTris->InsertNextCell(3);
                newTris->InsertCellPoint(centerPtId);
                newTris->InsertCellPoint(prevPtId);
                newTris->InsertCellPoint(firstPtId);

                // copy point data
                for( int i=0; i<numArrays; i++ )
                {
                    double data[9], data1[9], data2[9], d[9];
                    vtkDataArray* arr = inPD->GetArray(i);
                    arr->GetTuple(ptId1, data1);
                    arr->GetTuple(ptId2, data2);

                    for( int k=0; k<arr->GetNumberOfComponents(); k++) 
                        data[k] = 0.5*data1[k] + 0.5*data2[k];

                    for( int j=0; j<this->Samples+1; j++) 
                    {
                        outPD->GetArray(i)->InsertNextTuple(data);
                    }
                }
            }
        }

        this->UpdateProgress((float)cellId/numCells);
    }

    output->SetPoints(newPoints);
    output->SetPolys(newTris);
    newPoints->Delete();
    newTris->Delete();

    output->Squeeze();

    return 1;
}


//----------------------------------------------------------------------------
void ExtractBrokenBonds::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}


vtkSmartPointer<vtkIdList> ExtractBrokenBonds::getConnectedPoints( vtkPolyData * pd, vtkIdType ptId, vtkSmartPointer<vtkIdList> bondIds /*= NULL*/)
{
    vtkSmartPointer<vtkIdList> cellIds  = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> neighPts = vtkSmartPointer<vtkIdList>::New();

    pd->GetPointCells( ptId, cellIds);

    for(vtkIdType i = 0; i < cellIds->GetNumberOfIds(); i++)
    {
        vtkIdType cellId = cellIds->GetId(i);

        vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
        pd->GetCellPoints( cellId, pointIdList);

        vtkIdType ptIdN = pointIdList->GetId(0);
        if( ptIdN == ptId)
            ptIdN = pointIdList->GetId(1);

        neighPts->InsertNextId(ptIdN);

        if( bondIds )
            bondIds->InsertNextId(cellId);
    }

    return neighPts;
}