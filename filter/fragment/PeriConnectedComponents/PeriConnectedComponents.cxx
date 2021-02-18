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


#include "PeriConnectedComponents.h"

#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkMath.h>
#include <vtkIdList.h>

#include <vector>

using namespace std;

vtkStandardNewMacro(PeriConnectedComponents);

PeriConnectedComponents::PeriConnectedComponents()
{
    this->MaxDamage = 0.5;
    this->MaxBondlength = 0.1;

    // by default, process active point scalars
    this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
        vtkDataSetAttributes::SCALARS);

    // by default, process active cell scalars
    this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS,
        vtkDataSetAttributes::SCALARS);
}

PeriConnectedComponents::~PeriConnectedComponents()
{
}
   

int PeriConnectedComponents::RequestData(vtkInformation *vtkNotUsed(request),
	vtkInformationVector **inputVector,
	vtkInformationVector *outputVector)
{

	//this->DebugOn();

    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    
    // get the input and output
    vtkPolyData *input = vtkPolyData::SafeDownCast(
        inInfo->Get(vtkPolyData::DATA_OBJECT()));
    vtkPolyData *output = vtkPolyData::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkDataArray *dmg = this->GetInputArrayToProcess(0, inputVector);

    vtkFloatArray *refBondLengths = vtkFloatArray::SafeDownCast(
        this->GetInputArrayToProcess(1, inputVector));

    if (dmg == NULL || refBondLengths == NULL)
    {
        vtkErrorMacro("No input array.");
        return 0;
    }
    if (dmg->GetNumberOfComponents() == 0 || refBondLengths->GetNumberOfComponents() == 0 )
    {
        vtkErrorMacro("Input array must have at least one component.");
        return 0;
    }

    // estimate connected components
    input->BuildLinks();
    vector<bool> visited;
    int numPoints = input->GetNumberOfPoints();
    visited.resize(numPoints);

    vtkSmartPointer<vtkIntArray> labels = vtkSmartPointer<vtkIntArray>::New();
    labels->SetName("labels");
    labels->SetNumberOfComponents(1);
    labels->SetNumberOfValues(numPoints);

    for( int i=0; i<numPoints; i++)
    {
        visited[i]=false;
        labels->SetValue(i, 0);
    }
                
    // find first node with damage value lower than threshold
    int startNode = 0;
    while( startNode < visited.size() && dmg->GetTuple1(startNode) > this->MaxDamage ) {
        visited[startNode] = true;
        startNode++;
    }

    // breadth first search for connected components
    int currentLabel = 1;
    vector<vtkIdType> stack;
    double radius2 = MaxBondlength*MaxBondlength;

    while( startNode < visited.size() )
    {
        this->UpdateProgress (startNode/(float)(numPoints));

        stack.push_back(startNode);
        visited[startNode] = true;
        labels->SetValue( startNode, currentLabel);

        while( !stack.empty())
        {
            int currentNode = stack.back();
            stack.pop_back();
            double p1[3], p2[3];
            input->GetPoint(currentNode, p1);

            vtkSmartPointer<vtkIdList> bondsIds = vtkSmartPointer<vtkIdList>::New();
            vtkSmartPointer<vtkIdList> conn = getConnectedPoints(input, currentNode, bondsIds);

            for( int i = 0; i< conn->GetNumberOfIds(); i++)
            {
                int nodeId = conn->GetId(i);
                int bondId = bondsIds->GetId(i);

                // estimate length of bond
                //input->GetPoint(nodeId, p2);
                //double dist2 = vtkMath::Distance2BetweenPoints(p1,p2);
                float refBondLength = refBondLengths->GetValue(bondId);

                // check if node belongs to current fragment
                if( !visited[nodeId] && 
                     //dist2<radius2 && 
                     refBondLength < this->MaxBondlength &&
                     dmg->GetTuple1(nodeId) < this->MaxDamage)
                {
                    stack.push_back(nodeId);
                    labels->SetValue(nodeId, currentLabel);
                    visited[nodeId]=true;
                }
            }
        }
        // find next start node with damage value lower than threshold
        while(  startNode < visited.size() && 
            ( visited[startNode] || dmg->GetTuple1(startNode) > this->MaxDamage ))
        {
            visited[startNode] = true;
            startNode++;
        }

        currentLabel++;
    }

    output->ShallowCopy(input);
    output->GetPointData()->AddArray(labels);
    output->Squeeze();

	return 1;
}

//---------------------------------------------------------------------------
vtkSmartPointer<vtkIdList> PeriConnectedComponents::
    getConnectedPoints( vtkPolyData * pd, vtkIdType ptId, vtkSmartPointer<vtkIdList> bondIds /*= NULL*/)
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

//----------------------------------------------------------------------------
int PeriConnectedComponents::FillInputPortInformation(int, vtkInformation *info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
}


//==============================================================================
void PeriConnectedComponents::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}