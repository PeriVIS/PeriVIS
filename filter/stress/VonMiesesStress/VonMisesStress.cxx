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


#include "VonMisesStress.h"
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

vtkStandardNewMacro(VonMisesStress);

//==============================================================================
VonMisesStress::VonMisesStress()
{
  this->SetNumberOfInputPorts(1);
  
  // by default, process active point tensors
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                                    vtkDataSetAttributes::TENSORS);

}
//==============================================================================
VonMisesStress::~VonMisesStress()
{
  // clean up
}

//==============================================================================
int VonMisesStress::RequestData(
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

  vtkDataArray *inTensor;
  double tensor[9];
  vtkIdType numPts, inPtId, ptIncr, i;
  int j;

  
  vtkDebugMacro(<<"Calculating VonMises Stress");

  vtkPointData *outPD = output->GetPointData();
  inTensor = this->GetInputArrayToProcess(0, inputVector);
  numPts = input->GetNumberOfPoints();

  if ( !inTensor || numPts < 1 )
    {
    vtkErrorMacro(<<"No data!");
    return 1;
    }
  if (( inTensor->GetNumberOfComponents() != 4 ) &&
      ( inTensor->GetNumberOfComponents() != 9 ))
    {
    vtkErrorMacro("Input array must be a gradient tensor with 4 or 9 components.");
    return 0;
    }

  // allocate mem for output data arrays
  vtkSmartPointer<vtkDoubleArray> normArray = vtkSmartPointer<vtkDoubleArray>::New();
  normArray->SetName("Von Mises Stress");
  normArray->SetNumberOfComponents(1);
  normArray->SetNumberOfTuples(numPts);

  //
  // Traverse all Input points, calculate and store symmetric and antisymmetric tensors
  //

  double lambdaMax = 0.0;
  if( inTensor->GetNumberOfComponents() == 4 ){
  
    for (inPtId=0; inPtId < numPts; inPtId++)
      {

      inTensor->GetTuple(inPtId, tensor);

      double sig_x = tensor[0];
      double sig_y = tensor[3];
      double tau_xy = tensor[1];
      double sig_v = sqrt( sig_x*sig_x + sig_y*sig_y - sig_x*sig_y + 3*tau_xy*tau_xy );

      //copy tensors
      normArray->SetTupleValue(inPtId, &sig_v );
        
      }    
  } 
  else /* if( inTensor->GetNumberOfComponents() == 9 ) */
  {

    for (inPtId=0; inPtId < numPts; inPtId++)
      {

      inTensor->GetTuple(inPtId, tensor);

      double sig_x = tensor[0];
      double sig_y = tensor[4];
      double sig_z = tensor[8];
      double tau_xy = tensor[1];
      double tau_xz = tensor[2];
      double tau_yz = tensor[5];

      double sig_v = sqrt(   sig_x*sig_x
                           + sig_y*sig_y
                           + sig_z*sig_z
                           - sig_x*sig_y
                           - sig_x*sig_z
                           - sig_y*sig_z
                           + 3*( tau_xy*tau_xy 
                               + tau_xz*tau_xz 
                               + tau_yz*tau_yz )
                           );

      //copy tensors
      normArray->SetTupleValue(inPtId, &sig_v );

      }
  }

  vtkDebugMacro(<<"Calculated " << numPts <<" von Mises stress values");

  //
  // Update output and release memory
  //

  output->ShallowCopy(input);

  vtkPointData *pd = output->GetPointData();
  pd->AddArray(normArray);
 
  output->Squeeze();

  return 1;
}

//==============================================================================
void VonMisesStress::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
