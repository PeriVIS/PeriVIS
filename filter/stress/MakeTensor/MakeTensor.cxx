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


#include "MakeTensor.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkCell.h"
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
#include "vtkTransform.h"
#include "vtkSmartPointer.h"

//#include "linalg.h"

vtkStandardNewMacro(MakeTensor);

//==============================================================================
MakeTensor::MakeTensor()
{
  this->SetNumberOfInputPorts(1);
  
  // by default, process active point tensors
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                                    vtkDataSetAttributes::SCALARS);

}
//==============================================================================
MakeTensor::~MakeTensor()
{
  // clean up
}

//==============================================================================
int MakeTensor::RequestData(
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
  vtkDataSet *output = vtkDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDataArray *inSxx, *inSyy, *inSzz, *inSxy, *inSxz, *inSyz;
  double tensor[9];
  vtkIdType numPts, inPtId, ptIncr, i;
  int j;

  
  vtkDebugMacro(<<"Creating Tensor from input scalar arrays");

  vtkPointData *outPD = output->GetPointData();
  inSxx = this->GetInputArrayToProcess(0, inputVector);
  inSyy = this->GetInputArrayToProcess(1, inputVector);
  inSzz = this->GetInputArrayToProcess(2, inputVector);
  inSxy = this->GetInputArrayToProcess(3, inputVector);
  inSxz = this->GetInputArrayToProcess(4, inputVector);
  inSyz = this->GetInputArrayToProcess(5, inputVector);
  numPts = input->GetNumberOfPoints();

  if ( !inSxx || !inSyy || !inSzz || !inSxy || !inSxz || !inSyz || numPts < 1 )
    {
    vtkErrorMacro(<<"No data!");
    return 1;
    }

  // allocate mem for output data array
  vtkSmartPointer<vtkDoubleArray> tensorArray = vtkSmartPointer<vtkDoubleArray>::New();
  tensorArray->SetName("Tensor");
  tensorArray->SetNumberOfComponents(9);
  tensorArray->SetNumberOfTuples(numPts);

  //
  // Traverse all Input points, calculate and store tensors
  //

    for (inPtId=0; inPtId < numPts; inPtId++)
    {
        tensor[0] = inSxx->GetTuple(inPtId)[0];
        tensor[4] = inSyy->GetTuple(inPtId)[0];
        tensor[8] = inSzz->GetTuple(inPtId)[0];
        tensor[1] = tensor[3] = inSxy->GetTuple(inPtId)[0];
        tensor[2] = tensor[6] = inSxz->GetTuple(inPtId)[0];
        tensor[5] = tensor[7] = inSyz->GetTuple(inPtId)[0];
          
          //copy tensors
        tensorArray->SetTupleValue(inPtId, tensor );

      }

  vtkDebugMacro(<<"Calculated " << numPts <<" tensor norm values");

  //
  // Update output and release memory
  //

  output->ShallowCopy(input);

  vtkPointData *pd = output->GetPointData();
  pd->AddArray(tensorArray);
 
  output->Squeeze();

  return 1;
}

//==============================================================================
void MakeTensor::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
