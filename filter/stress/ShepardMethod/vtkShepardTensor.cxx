/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkShepardTensor.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkShepardTensor.h"

#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"

vtkStandardNewMacro(vtkShepardTensor);

// Construct with sample dimensions=(50,50,50) and so that model bounds are
// automatically computed from input. Null value for each unvisited output
// point is 0.0. Maximum distance is 0.25.
vtkShepardTensor::vtkShepardTensor()
{
  this->MaximumDistance = 0.25;

  this->ModelBounds[0] = -1.0;
  this->ModelBounds[1] =  1.0;
  this->ModelBounds[2] = -1.0;
  this->ModelBounds[3] =  1.0;
  this->ModelBounds[4] = -1.0;
  this->ModelBounds[5] =  1.0;
  this->UseInputBounds= true;
  this->ProcessAllDataArrays = false;

  this->SampleDimensions[0] = 50;
  this->SampleDimensions[1] = 50;
  this->SampleDimensions[2] = 50;

  this->NullValue = 0.0;

  this->Dimension = 3;
  
  this->SetNumberOfInputPorts(1);
  
  // by default, process active point tensors
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::TENSORS);
}

// Compute ModelBounds from input geometry.
double vtkShepardTensor::ComputeModelBounds(double origin[3],
                                            double spacing[3])
{
  double *bounds, maxDist;
  int i, adjustBounds=0;

  // compute model bounds if not set previously
  if ( this->UseInputBounds ||
	   this->ModelBounds[0] >= this->ModelBounds[1] ||
       this->ModelBounds[2] >= this->ModelBounds[3] ||
       this->ModelBounds[4] >= this->ModelBounds[5] )
    {
    adjustBounds = 1;
    vtkDataSet *ds = vtkDataSet::SafeDownCast(this->GetInput());
    // ds better be non null otherwise something is very wrong here
    bounds = ds->GetBounds();
    }
  else
    {
    bounds = this->ModelBounds;
    }

  for (maxDist=0.0, i=0; i<this->Dimension; i++)
    {
    if ( (bounds[2*i+1] - bounds[2*i]) > maxDist )
      {
      maxDist = bounds[2*i+1] - bounds[2*i];
      }
    }
  maxDist *= this->MaximumDistance;

  // adjust bounds so model fits strictly inside (only if not set previously)
  if ( adjustBounds )
  {
	if(this->UseInputBounds)
	{
		for (i=0; i<this->Dimension; i++)
		{
			this->ModelBounds[2*i] = bounds[2*i];
			this->ModelBounds[2*i+1] = bounds[2*i+1];
		}
	}
	else
	{
		for (i=0; i<this->Dimension; i++)
		{
		  this->ModelBounds[2*i] = bounds[2*i] - maxDist;
		  this->ModelBounds[2*i+1] = bounds[2*i+1] + maxDist;
		}
    }
  }

  // Set volume origin and data spacing
  for (i=0; i<this->Dimension; i++)
    {
    origin[i] = this->ModelBounds[2*i];
    spacing[i] = (this->ModelBounds[2*i+1] - this->ModelBounds[2*i])
            / (this->SampleDimensions[i] - 1);
    }
  if( this->Dimension == 2) {
      origin[2] = 0.0;
      spacing[i] = 0.0;
  }
  return maxDist;
}

int vtkShepardTensor::RequestInformation (
  vtkInformation * vtkNotUsed(request),
  vtkInformationVector ** vtkNotUsed( inputVector ),
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  int i;
  double ar[3]={0,0,0}, origin[3]={0,0,0};

  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
               0, this->SampleDimensions[0]-1,
               0, this->SampleDimensions[1]-1,
               0, this->SampleDimensions[2]-1);

  for (i=0; i < this->Dimension; i++)
    {
    origin[i] = this->ModelBounds[2*i];
    if ( this->SampleDimensions[i] <= 1 )
      {
      ar[i] = 1;
      }
    else
      {
      ar[i] = (this->ModelBounds[2*i+1] - this->ModelBounds[2*i])
              / (this->SampleDimensions[i] - 1);
      }
    }
  outInfo->Set(vtkDataObject::ORIGIN(),origin,3);
  outInfo->Set(vtkDataObject::SPACING(),ar,3);

  return 1;
}

int vtkShepardTensor::RequestData(
  vtkInformation* vtkNotUsed( request ),
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // get the input
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkDataSet *input = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  // get the output
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkImageData *output = vtkImageData::SafeDownCast( 
      outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkIdType ptId, i;
  int j, k;
  double *px, x[3], s[9], *sum, spacing[3], origin[3];

  double maxDistance, distance2, inTensor[9];
  vtkDataArray *inTensors;
  vtkIdType numPts, numNewPts, idx;
  int min[3], max[3];
  int jkFactor;

  vtkDebugMacro(<< "Executing Shepard method");

  // Check input
  //
  if ( (numPts=input->GetNumberOfPoints()) < 1 )
    {
    vtkErrorMacro(<<"Points must be defined!");
    return 1;
    }
  
  inTensors = this->GetInputArrayToProcess(0, inputVector);
  
  if (inTensors == NULL )
    {
    vtkErrorMacro(<<"Tensor must be defined!");
    return 1;
    }

  int numComponents = inTensors->GetNumberOfComponents();

  // Allocate
  //
  numNewPts = this->SampleDimensions[0] * this->SampleDimensions[1] * this->SampleDimensions[2];
  
  vtkFloatArray* newTensors = vtkFloatArray::New();
  newTensors->SetNumberOfComponents(numComponents);
  newTensors->SetNumberOfTuples(numNewPts);

  vtkPointData *inPD = input->GetPointData(),
               *outPD = output->GetPointData();

  // copy point data arrays
  int numArrays = inPD->GetNumberOfArrays();

  if( ProcessAllDataArrays ) {
      for( int i=0; i<numArrays; i++ ){
          int type = inPD->GetArray(i)->GetDataType();
          vtkDataArray *arr = vtkDataArray::CreateDataArray(type);
          arr->SetNumberOfComponents(inPD->GetArray(i)->GetNumberOfComponents());
          arr->SetName(inPD->GetArray(i)->GetName());
          arr->SetNumberOfTuples(numNewPts);
          outPD->AddArray(arr);
          arr->Delete();
      }
  }

  newTensors->SetName( inTensors->GetName());

  const double nullTensor[9] = { 0.0, 0.0, 0.0, 
                                 0.0, 0.0, 0.0, 
                                 0.0, 0.0, 0.0 };
  const double maxTensor[9] = { VTK_FLOAT_MAX, VTK_FLOAT_MAX, VTK_FLOAT_MAX, 
                                VTK_FLOAT_MAX, VTK_FLOAT_MAX, VTK_FLOAT_MAX,
                                VTK_FLOAT_MAX, VTK_FLOAT_MAX, VTK_FLOAT_MAX };

  sum = new double[numNewPts];
  for (i=0; i<numNewPts; i++)
    {
        if( ProcessAllDataArrays ) {
            for( int ai=0; ai<numArrays; ai++ ){
                outPD->GetArray(ai)->SetTuple(i, nullTensor);
            }
        }
        else {
            newTensors->SetTuple(i, nullTensor);
        }
        sum[i] = 0.0;
    }

  maxDistance = this->ComputeModelBounds(origin,spacing);


  // Traverse all input points.
  // Each input point affects voxels within maxDistance.
  //
  for (ptId=0; ptId < numPts; ptId++)
    {
    if ( ! (ptId % 1000) )
      {
      vtkDebugMacro(<<"Inserting point #" << ptId);
      this->UpdateProgress (ptId/(double)numPts);
      if (this->GetAbortExecute())
        {
        break;
        }
      }

    px = input->GetPoint(ptId);
    inTensors->GetTuple(ptId, inTensor);

    for (i=0; i<this->Dimension; i++) //compute dimensional bounds in data set
      {
      double amin = static_cast<double>(
        (px[i] - maxDistance) - origin[i]) / spacing[i];
      double amax = static_cast<double>(
        (px[i] + maxDistance) - origin[i]) / spacing[i];
      min[i] = static_cast<int>(amin);
      max[i] = static_cast<int>(amax);

      if (min[i] < amin)
        {
        min[i]++; // round upward to nearest integer to get min[i]
        }
      if (max[i] > amax)
        {
        max[i]--; // round downward to nearest integer to get max[i]
        }

      if (min[i] < 0)
        {
        min[i] = 0; // valid range check
        }
      if (max[i] >= this->SampleDimensions[i])
        {
        max[i] = this->SampleDimensions[i] - 1;
        }
      }

    for (i=0; i<=this->Dimension; i++) //compute dimensional bounds in data set
      {
      min[i] = static_cast<int>(
        static_cast<double>((px[i] - maxDistance) - origin[i]) / spacing[i]);
      max[i] = static_cast<int>(
        static_cast<double>((px[i] + maxDistance) - origin[i]) / spacing[i]);
      if (min[i] < 0)
        {
        min[i] = 0;
        }
      if (max[i] >= this->SampleDimensions[i])
        {
        max[i] = this->SampleDimensions[i] - 1;
        }
      }

    if( this->Dimension == 2 )
      {
          //jkFactor = this->SampleDimensions[0]*this->SampleDimensions[1];
          x[2] = 0.0;
          for (j = min[1]; j <= max[1]; j++)
          {
              x[1] = spacing[1] * j + origin[1];
              for (i = min[0]; i <= max[0]; i++)
              {
                  x[0] = spacing[0] * i + origin[0];
                  idx = this->SampleDimensions[0]*j + i;

                  distance2 = vtkMath::Distance2BetweenPoints(x,px);

                  if ( distance2 == 0.0 )
                  {
                      sum[idx] = VTK_DOUBLE_MAX;
                      if( ProcessAllDataArrays ) {
                          for( int ai=0; ai<numArrays; ai++ ){
                              outPD->GetArray(ai)->SetTuple(idx, maxTensor);
                          }
                      } else {
                          newTensors->SetTuple(idx, maxTensor);
                      }
                  }
                  else
                  {
                      sum[idx] += 1.0 / distance2;

                      if( ProcessAllDataArrays ) {
                          for( int ai=0; ai<numArrays; ai++ ){
                              outPD->GetArray(ai)->GetTuple(idx, s);
                              inPD->GetArray(ai)->GetTuple(ptId, inTensor);
                              for( int i_s = 0; i_s<outPD->GetArray(ai)->GetNumberOfComponents(); i_s++ )
                              {
                                  s[i_s] = s[i_s]+(inTensor[i_s]/distance2);
                              }
                              outPD->GetArray(ai)->SetTuple(idx, s);
                          }
                      } else {
                          newTensors->GetTuple(idx, s);
                          for( int i_s = 0; i_s<numComponents; i_s++ )
                          {
                              //s = s+(inScalar/distance2)
                              s[i_s] = s[i_s]+(inTensor[i_s]/distance2);
                          }
                          newTensors->SetTuple(idx, s);
                      }
                  }
              }
          }
      }
    else /*if( this->Dimension == 3 )*/
      {
          jkFactor = this->SampleDimensions[0]*this->SampleDimensions[1];
          for (k = min[2]; k <= max[2]; k++)
          {
              x[2] = spacing[2] * k + origin[2];
              for (j = min[1]; j <= max[1]; j++)
              {
                  x[1] = spacing[1] * j + origin[1];
                  for (i = min[0]; i <= max[0]; i++)
                  {
                      x[0] = spacing[0] * i + origin[0];
                      idx = jkFactor*k + this->SampleDimensions[0]*j + i;

                      distance2 = vtkMath::Distance2BetweenPoints(x,px);

                      if ( distance2 == 0.0 )
                      {
                          sum[idx] = VTK_DOUBLE_MAX;
                          if( ProcessAllDataArrays) {
                              for( int ai=0; ai<numArrays; ai++ ){
                                  outPD->GetArray(ai)->SetTuple(idx, maxTensor);
                              }
                          } else {
                              newTensors->SetTuple(idx, maxTensor);
                          }
                      }
                      else
                      {
                          sum[idx] += 1.0 / distance2;

                          if( ProcessAllDataArrays ) {
                              for( int ai=0; ai<numArrays; ai++ ){
                                  outPD->GetArray(ai)->GetTuple(idx, s);
                                  inPD->GetArray(ai)->GetTuple(ptId, inTensor);
                                  for( int i_s = 0; i_s<outPD->GetArray(ai)->GetNumberOfComponents(); i_s++ )
                                  {
                                      s[i_s] = s[i_s]+(inTensor[i_s]/distance2);
                                  }
                                  outPD->GetArray(ai)->SetTuple(idx, s);
                              }
                          } else {
                              newTensors->GetTuple(idx, s);
                              for( int i_s = 0; i_s<numComponents; i_s++ )
                              {
                                  //s = s+(inScalar/distance2)
                                  s[i_s] = s[i_s]+(inTensor[i_s]/distance2);
                              }
                              newTensors->SetTuple(idx, s);
                          }
                      }
                  }
              }
          }   
      }


    }

  // Run through scalars and compute final values
  //
  for (ptId=0; ptId<numNewPts; ptId++)
    {
    newTensors->GetTuple(ptId, s);
    if ( sum[ptId] != 0.0 )
      {
          if( ProcessAllDataArrays ) {
              for( int ai=0; ai<numArrays; ai++ )
              {
                  outPD->GetArray(ai)->GetTuple(ptId, s);
                  for( int i_s = 0; i_s<outPD->GetArray(ai)->GetNumberOfComponents(); i_s++ )
                  {
                      s[i_s] = s[i_s] / sum[ptId];
                  }
                  outPD->GetArray(ai)->SetTuple(ptId, s);
              }
          } else {
              for( int idx = 0; idx<numComponents; idx++ ) 
                {
                s[idx] = s[idx] / sum[ptId];
                }
              newTensors->SetTuple(ptId, s);
          }
      }
    else
      {
          if( ProcessAllDataArrays) {
              for( int ai=0; ai<numArrays; ai++ ){
                  outPD->GetArray(ai)->SetTuple(idx, nullTensor);
              }
          } else {
              newTensors->SetTuple(ptId,nullTensor);
          }
      }
    }

  output->SetExtent( 0, this->SampleDimensions[0]-1,
                     0, this->SampleDimensions[1]-1,
                     0, this->SampleDimensions[2]-1 );
  output->SetOrigin(origin);
  output->SetSpacing(spacing);
  
  if( !this->ProcessAllDataArrays) {
    output->GetPointData()->AddArray(newTensors);
  }

  newTensors->Delete();

  // Update self
  //
  delete [] sum;
  return 1;
}

// Set the i-j-k dimensions on which to sample the distance function.
void vtkShepardTensor::SetSampleDimensions(int i, int j, int k)
{
  int dim[3];

  dim[0] = i;
  dim[1] = j;
  dim[2] = k;

  this->SetSampleDimensions(dim);
}

// Set the i-j-k dimensions on which to sample the distance function.
void vtkShepardTensor::SetSampleDimensions(int dim[3])
{
  int dataDim, i;

  vtkDebugMacro(<< " setting SampleDimensions to (" << dim[0] << ","
                << dim[1] << "," << dim[2] << ")");

  if ( dim[0] != this->SampleDimensions[0] ||
       dim[1] != this->SampleDimensions[1] ||
       dim[2] != this->SampleDimensions[2] )
    {
    if ( dim[0]<1 || dim[1]<1 || dim[2]<1 )
      {
      vtkErrorMacro (<< "Bad Sample Dimensions, retaining previous values");
      return;
      }

    for (dataDim=0, i=0; i<3 ; i++)
      {
      if (dim[i] > 1)
        {
        dataDim++;
        }
      }

    this->Dimension = dataDim;

    if ( dataDim  < 2 )
      {
      vtkErrorMacro(<<"Sample dimensions must define a volume!");
      return;
      }

    for ( i=0; i<3; i++)
      {
      this->SampleDimensions[i] = dim[i];
      }

    this->Modified();
    }
}

int vtkShepardTensor::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}

void vtkShepardTensor::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Maximum Distance: " << this->MaximumDistance << "\n";

  os << indent << "Sample Dimensions: (" << this->SampleDimensions[0] << ", "
               << this->SampleDimensions[1] << ", "
               << this->SampleDimensions[2] << ")\n";

  os << indent << "ModelBounds: \n";
  os << indent << "  Xmin,Xmax: (" << this->ModelBounds[0] << ", " << this->ModelBounds[1] << ")\n";
  os << indent << "  Ymin,Ymax: (" << this->ModelBounds[2] << ", " << this->ModelBounds[3] << ")\n";
  os << indent << "  Zmin,Zmax: (" << this->ModelBounds[4] << ", " << this->ModelBounds[5] << ")\n";

  os << indent << "Null Value: " << this->NullValue << "\n";

}
