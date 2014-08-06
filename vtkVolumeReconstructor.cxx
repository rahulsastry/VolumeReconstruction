/*=Plus=header=begin======================================================
  Program: Plus
  Copyright (c) Laboratory for Percutaneous Surgery. All rights reserved.
  See License.txt for details.
=========================================================Plus=header=end*/

#include "PlusConfigure.h"

#include "vtkVolumeReconstructor.h"

#include <limits>

#include "vtkImageImport.h" 
#include "vtkImageData.h" 
#include "vtkImageViewer.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
#include "vtkTransform.h"
#include "vtkXMLUtilities.h"
#include "vtkImageExtractComponents.h"
#include "vtkDataSetWriter.h"

#include "vtkPasteSliceIntoVolume.h"
#include "vtkFillHolesInVolume.h"
#include "vtkTrackedFrameList.h"
#include "TrackedFrame.h"
#include "vtkTransformRepository.h"

#include "metaImage.h"

#include "vtkMath.h"
#include <cmath>
#include <vector>
#include <iostream>

vtkCxxRevisionMacro(vtkVolumeReconstructor, "$Revisions: 1.0 $");
vtkStandardNewMacro(vtkVolumeReconstructor);

//----------------------------------------------------------------------------
vtkVolumeReconstructor::vtkVolumeReconstructor()
: ImageCoordinateFrame(NULL)
, ReferenceCoordinateFrame(NULL)
{
  this->ReconstructedVolume = vtkSmartPointer<vtkImageData>::New();
  this->Reconstructor = vtkPasteSliceIntoVolume::New();  
  this->HoleFiller = vtkFillHolesInVolume::New();  
  this->FillHoles = false;
  this->SkipInterval = 1;
  this->ReconstructedVolumeUpdatedTime = 0;
}

//----------------------------------------------------------------------------
vtkVolumeReconstructor::~vtkVolumeReconstructor()
{
  if (this->Reconstructor)
  {
    this->Reconstructor->Delete();
    this->Reconstructor=NULL;
  }

  if (this->HoleFiller)
  {
    this->HoleFiller->Delete();
    this->HoleFiller=NULL;
  }
  SetImageCoordinateFrame(NULL);
  SetReferenceCoordinateFrame(NULL);
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
PlusStatus vtkVolumeReconstructor::ReadConfiguration(vtkXMLDataElement* config)
{
  XML_FIND_NESTED_ELEMENT_REQUIRED(reconConfig, config, "VolumeReconstruction");

  XML_READ_STRING_ATTRIBUTE_OPTIONAL(ReferenceCoordinateFrame, reconConfig);
  XML_READ_STRING_ATTRIBUTE_OPTIONAL(ImageCoordinateFrame, reconConfig);

  XML_READ_VECTOR_ATTRIBUTE_REQUIRED(double, 3, OutputSpacing, reconConfig);
  XML_READ_VECTOR_ATTRIBUTE_OPTIONAL(double, 3, OutputOrigin, reconConfig);
  XML_READ_VECTOR_ATTRIBUTE_OPTIONAL(int, 6, OutputExtent, reconConfig);

  XML_READ_VECTOR_ATTRIBUTE_OPTIONAL(int, 2, ClipRectangleOrigin, reconConfig);
  XML_READ_VECTOR_ATTRIBUTE_OPTIONAL(int, 2, ClipRectangleSize, reconConfig);

  XML_READ_VECTOR_ATTRIBUTE_OPTIONAL(double, 2, FanAngles, reconConfig);  // DEPRECATED (replaced by FanAnglesDeg)
  XML_READ_VECTOR_ATTRIBUTE_OPTIONAL(double, 2, FanAnglesDeg, reconConfig);
  
  XML_READ_VECTOR_ATTRIBUTE_OPTIONAL(double, 2, FanOrigin, reconConfig);  // DEPRECATED (replaced by FanOriginPixel)
  XML_READ_VECTOR_ATTRIBUTE_OPTIONAL(double, 2, FanOriginPixel, reconConfig);

  XML_READ_SCALAR_ATTRIBUTE_OPTIONAL(double, FanDepth, reconConfig);   // DEPRECATED (replaced by FanRadiusStopPixel)
  XML_READ_SCALAR_ATTRIBUTE_OPTIONAL(double, FanRadiusStartPixel, reconConfig);
  XML_READ_SCALAR_ATTRIBUTE_OPTIONAL(double, FanRadiusStopPixel, reconConfig);

  XML_READ_SCALAR_ATTRIBUTE_OPTIONAL(int, SkipInterval, reconConfig);
  if (this->SkipInterval < 1)
  {
    LOG_WARNING("SkipInterval in the config file must be greater or equal to 1. Resetting to 1");
    SkipInterval = 1;
  }

  // reconstruction options
  XML_READ_ENUM2_ATTRIBUTE_OPTIONAL(Interpolation, reconConfig, \
    this->Reconstructor->GetInterpolationModeAsString(vtkPasteSliceIntoVolume::LINEAR_INTERPOLATION), vtkPasteSliceIntoVolume::LINEAR_INTERPOLATION, \
    this->Reconstructor->GetInterpolationModeAsString(vtkPasteSliceIntoVolume::NEAREST_NEIGHBOR_INTERPOLATION), vtkPasteSliceIntoVolume::NEAREST_NEIGHBOR_INTERPOLATION);

  XML_READ_ENUM2_ATTRIBUTE_OPTIONAL(Calculation, reconConfig, \
    this->Reconstructor->GetCalculationModeAsString(vtkPasteSliceIntoVolume::WEIGHTED_AVERAGE), vtkPasteSliceIntoVolume::WEIGHTED_AVERAGE, \
    this->Reconstructor->GetCalculationModeAsString(vtkPasteSliceIntoVolume::MAXIMUM), vtkPasteSliceIntoVolume::MAXIMUM);

  XML_READ_ENUM3_ATTRIBUTE_OPTIONAL(Optimization, reconConfig, \
    this->Reconstructor->GetOptimizationModeAsString(vtkPasteSliceIntoVolume::FULL_OPTIMIZATION), vtkPasteSliceIntoVolume::FULL_OPTIMIZATION, \
    this->Reconstructor->GetOptimizationModeAsString(vtkPasteSliceIntoVolume::PARTIAL_OPTIMIZATION), vtkPasteSliceIntoVolume::PARTIAL_OPTIMIZATION, \
    this->Reconstructor->GetOptimizationModeAsString(vtkPasteSliceIntoVolume::NO_OPTIMIZATION), vtkPasteSliceIntoVolume::NO_OPTIMIZATION);

  XML_READ_ENUM2_ATTRIBUTE_OPTIONAL(Compounding, reconConfig, "ON", true, "OFF", false);

  XML_READ_SCALAR_ATTRIBUTE_OPTIONAL(int, NumberOfThreads, reconConfig);

  XML_READ_ENUM2_ATTRIBUTE_OPTIONAL(FillHoles, reconConfig, "ON", true, "OFF", false);

  // Find and read kernels. First for loop counts the number of kernels to allocate, second for loop stores them
  if (this->FillHoles) 
  {
    // load input for kernel size, stdev, etc...
    XML_FIND_NESTED_ELEMENT_REQUIRED(holeFilling, reconConfig, "HoleFilling");
    if (this->HoleFiller->ReadConfiguration(holeFilling)!=PLUS_SUCCESS)
    {
      return PLUS_FAIL;
    }
  }

  this->Modified();

  return PLUS_SUCCESS;
}

//----------------------------------------------------------------------------
// Get the XML element describing the freehand object
PlusStatus vtkVolumeReconstructor::WriteConfiguration(vtkXMLDataElement *config)
{
  XML_FIND_NESTED_ELEMENT_CREATE_IF_MISSING(reconConfig, config, "VolumeReconstruction");

  reconConfig->SetAttribute("ImageCoordinateFrame", this->ImageCoordinateFrame);
  reconConfig->SetAttribute("ReferenceCoordinateFrame", this->ReferenceCoordinateFrame);

  // output parameters
  reconConfig->SetVectorAttribute("OutputSpacing", 3, this->Reconstructor->GetOutputSpacing());
  reconConfig->SetVectorAttribute("OutputOrigin", 3, this->Reconstructor->GetOutputOrigin());
  reconConfig->SetVectorAttribute("OutputExtent", 6, this->Reconstructor->GetOutputExtent());

  // clipping parameters
  reconConfig->SetVectorAttribute("ClipRectangleOrigin", 2, this->Reconstructor->GetClipRectangleOrigin());
  reconConfig->SetVectorAttribute("ClipRectangleSize", 2, this->Reconstructor->GetClipRectangleSize());

  // Fan parameters
  // remove deprecated attributes
  XML_REMOVE_ATTRIBUTE(reconConfig, "FanDepth");
  XML_REMOVE_ATTRIBUTE(reconConfig, "FanOrigin");
  XML_REMOVE_ATTRIBUTE(reconConfig, "FanAngles");
  if (this->Reconstructor->FanClippingApplied())
  {
    reconConfig->SetVectorAttribute("FanAnglesDeg", 2, this->Reconstructor->GetFanAnglesDeg());
    // Image spacing is 1.0, so reconstructor's fan origin and radius values are in pixels
    reconConfig->SetVectorAttribute("FanOriginPixel", 2, this->Reconstructor->GetFanOrigin());
    reconConfig->SetDoubleAttribute("FanRadiusStartPixel", this->Reconstructor->GetFanRadiusStart());
    reconConfig->SetDoubleAttribute("FanRadiusStopPixel", this->Reconstructor->GetFanRadiusStop());
  }
  else
  {
    XML_REMOVE_ATTRIBUTE(reconConfig, "FanAnglesDeg");
    XML_REMOVE_ATTRIBUTE(reconConfig, "FanOriginPixel");
    XML_REMOVE_ATTRIBUTE(reconConfig, "FanRadiusStartPixel");
    XML_REMOVE_ATTRIBUTE(reconConfig, "FanRadiusStopPixel");
  }

  // reconstruction options
  reconConfig->SetAttribute("Interpolation", this->Reconstructor->GetInterpolationModeAsString(this->Reconstructor->GetInterpolationMode()));
  reconConfig->SetAttribute("Optimization", this->Reconstructor->GetOptimizationModeAsString(this->Reconstructor->GetOptimization()));
  reconConfig->SetAttribute("Compounding", this->Reconstructor->GetCompounding()?"On":"Off");

  if (this->Reconstructor->GetNumberOfThreads()>0)
  {
    reconConfig->SetIntAttribute("NumberOfThreads", this->Reconstructor->GetNumberOfThreads());
  }
  else
  {
    XML_REMOVE_ATTRIBUTE(reconConfig, "NumberOfThreads");
  }

  return PLUS_SUCCESS;
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::AddImageToExtent( vtkImageData *image, vtkMatrix4x4* imageToReference, double* extent_Ref)
{
  // Output volume is in the Reference coordinate system.

  // Prepare the four corner points of the input US image.
  int* frameExtent=image->GetExtent();
  std::vector< double* > corners_ImagePix;
  double minX=frameExtent[0];
  double maxX=frameExtent[1];
  double minY=frameExtent[2];
  double maxY=frameExtent[3];
  if (this->Reconstructor->GetClipRectangleSize()[0]>0 && this->Reconstructor->GetClipRectangleSize()[1]>0)
  {
    // Clipping rectangle is specified
    minX=std::max<double>(minX, this->Reconstructor->GetClipRectangleOrigin()[0]);
    maxX=std::min<double>(maxX, this->Reconstructor->GetClipRectangleOrigin()[0]+this->Reconstructor->GetClipRectangleSize()[0]);
    minY=std::max<double>(minY, this->Reconstructor->GetClipRectangleOrigin()[1]);
    maxY=std::min<double>(maxY, this->Reconstructor->GetClipRectangleOrigin()[1]+this->Reconstructor->GetClipRectangleSize()[1]);
  }
  double c0[ 4 ] = { minX, minY, 0,  1 };
  double c1[ 4 ] = { minX, maxY, 0,  1 };
  double c2[ 4 ] = { maxX, minY, 0,  1 };
  double c3[ 4 ] = { maxX, maxY, 0,  1 };
  corners_ImagePix.push_back( c0 );
  corners_ImagePix.push_back( c1 );
  corners_ImagePix.push_back( c2 );
  corners_ImagePix.push_back( c3 );

  // Transform the corners to Reference and expand the extent if needed
  for ( unsigned int corner = 0; corner < corners_ImagePix.size(); ++corner )
  {
    double corner_Ref[ 4 ] = { 0, 0, 0, 1 }; // position of the corner in the Reference coordinate system
    imageToReference->MultiplyPoint( corners_ImagePix[corner], corner_Ref );

    for ( int axis = 0; axis < 3; axis ++ )
    {
      if ( corner_Ref[axis] < extent_Ref[axis*2] )
      {
        // min extent along this coord axis has to be decreased
        extent_Ref[axis*2]=corner_Ref[axis];
      }
      if ( corner_Ref[axis] > extent_Ref[axis*2+1] )
      {
        // max extent along this coord axis has to be increased
        extent_Ref[axis*2+1]=corner_Ref[axis];
      }
    }
  }
} 

//----------------------------------------------------------------------------
PlusStatus vtkVolumeReconstructor::GetImageToReferenceTransformName(PlusTransformName& imageToReferenceTransformName)
{
  if (this->ReferenceCoordinateFrame!=NULL && this->ImageCoordinateFrame!=NULL)
  {
    // image to reference transform is specified in the XML tree
    imageToReferenceTransformName=PlusTransformName(this->ImageCoordinateFrame,this->ReferenceCoordinateFrame);
    if (!imageToReferenceTransformName.IsValid())
    { 
      LOG_ERROR("Failed to set ImageToReference transform name from '" << this->ImageCoordinateFrame <<"' to '"<<this->ReferenceCoordinateFrame<<"'" ); 
      return PLUS_FAIL;
    }
    return PLUS_SUCCESS;
  }
  if (this->ImageCoordinateFrame==NULL)
  {
    LOG_ERROR("Image coordinate frame name is undefined");
  }
  if (this->ReferenceCoordinateFrame==NULL)
  {
    LOG_ERROR("Reference coordinate frame name is undefined");
  }
  return PLUS_FAIL;
}

//----------------------------------------------------------------------------
PlusStatus vtkVolumeReconstructor::SetOutputExtentFromFrameList(vtkTrackedFrameList* trackedFrameList, vtkTransformRepository* transformRepository)
{
  PlusTransformName imageToReferenceTransformName;
  if (GetImageToReferenceTransformName(imageToReferenceTransformName)!=PLUS_SUCCESS)
  {
    LOG_ERROR("Invalid ImageToReference transform name"); 
    return PLUS_FAIL; 
  }

  if ( trackedFrameList == NULL )
  {
    LOG_ERROR("Failed to set output extent from tracked frame list - input frame list is NULL!"); 
    return PLUS_FAIL; 
  }
  if ( trackedFrameList->GetNumberOfTrackedFrames() == 0)
  {
    LOG_ERROR("Failed to set output extent from tracked frame list - input frame list is empty!"); 
    return PLUS_FAIL; 
  }

  if ( transformRepository == NULL )
  {
    LOG_ERROR("Failed to set output extent from tracked frame list - input transform repository is NULL!"); 
    return PLUS_FAIL; 
  }

  double extent_Ref[6]=
  {
    VTK_DOUBLE_MAX, VTK_DOUBLE_MIN,
    VTK_DOUBLE_MAX, VTK_DOUBLE_MIN,
    VTK_DOUBLE_MAX, VTK_DOUBLE_MIN
  };

  const int numberOfFrames = trackedFrameList->GetNumberOfTrackedFrames();
  int numberOfValidFrames = 0;
  for (int frameIndex = 0; frameIndex < numberOfFrames; ++frameIndex )
  {
    TrackedFrame* frame = trackedFrameList->GetTrackedFrame( frameIndex );
    
    if ( transformRepository->SetTransforms(*frame) != PLUS_SUCCESS )
    {
      LOG_ERROR("Failed to update transform repository with tracked frame!"); 
      return PLUS_FAIL; 
    }

    // Get transform
    bool isMatrixValid(false); 
    vtkSmartPointer<vtkMatrix4x4> imageToReferenceTransformMatrix=vtkSmartPointer<vtkMatrix4x4>::New();
    if ( transformRepository->GetTransform(imageToReferenceTransformName, imageToReferenceTransformMatrix, &isMatrixValid ) != PLUS_SUCCESS )
    {
      std::string strImageToReferenceTransformName; 
      imageToReferenceTransformName.GetTransformName(strImageToReferenceTransformName); 
      LOG_ERROR("Failed to get transform '"<<strImageToReferenceTransformName<<"' from transform repository!"); 
      return PLUS_FAIL; 
    }

    if ( isMatrixValid )
    {
      numberOfValidFrames++;

      // Get image (only the frame extents will be used)
      vtkImageData* frameImage=trackedFrameList->GetTrackedFrame(frameIndex)->GetImageData()->GetImage();

      // Expand the extent_Ref to include this frame
      AddImageToExtent(frameImage, imageToReferenceTransformMatrix, extent_Ref);
    }
  }

  LOG_DEBUG("Automatic volume extent computation from frames used "<<numberOfValidFrames<<" out of "<<numberOfFrames<<" (probably wrong image or reference coordinate system was defined or all transforms were invalid)");
  if (numberOfValidFrames==0)
  {
    std::string strImageToReferenceTransformName; 
    imageToReferenceTransformName.GetTransformName(strImageToReferenceTransformName);
    LOG_ERROR("Automatic volume extent computation failed, there were no valid "<<strImageToReferenceTransformName<<" transform available in the whole sequence");
    return PLUS_FAIL;
  }

  // Set the output extent from the current min and max values, using the user-defined image resolution.
  int outputExtent[ 6 ] = { 0, 0, 0, 0, 0, 0 };
  double* outputSpacing = this->Reconstructor->GetOutputSpacing();
  outputExtent[ 1 ] = int( ( extent_Ref[1] - extent_Ref[0] ) / outputSpacing[ 0 ] );
  outputExtent[ 3 ] = int( ( extent_Ref[3] - extent_Ref[2] ) / outputSpacing[ 1 ] );
  outputExtent[ 5 ] = int( ( extent_Ref[5] - extent_Ref[4] ) / outputSpacing[ 2 ] );

  this->Reconstructor->SetOutputScalarMode(trackedFrameList->GetTrackedFrame(0)->GetImageData()->GetImage()->GetScalarType());
  this->Reconstructor->SetOutputExtent( outputExtent );
  this->Reconstructor->SetOutputOrigin( extent_Ref[0], extent_Ref[2], extent_Ref[4] ); 
  try
  {
    if (this->Reconstructor->ResetOutput()!=PLUS_SUCCESS) // :TODO: call this automatically
    {
      LOG_ERROR("Failed to initialize output of the reconstructor");
      return PLUS_FAIL;
    }
  }
  catch(std::bad_alloc& e)
  {
    cerr << e.what() << endl;
    LOG_ERROR("StartReconstruction failed due to out of memory. Try to reduce the size or increase spacing of the output volume.");
    return PLUS_FAIL;
  }

  this->Modified();

  return PLUS_SUCCESS;
}

//----------------------------------------------------------------------------
PlusStatus vtkVolumeReconstructor::AddTrackedFrame(TrackedFrame* frame, vtkTransformRepository* transformRepository, bool* insertedIntoVolume/*=NULL*/)
{
  PlusTransformName imageToReferenceTransformName;
  if (GetImageToReferenceTransformName(imageToReferenceTransformName)!=PLUS_SUCCESS)
  {
    LOG_ERROR("Invalid ImageToReference transform name"); 
    return PLUS_FAIL; 
  }

  if ( frame == NULL )
  {
    LOG_ERROR("Failed to add tracked frame to volume - input frame is NULL"); 
    return PLUS_FAIL; 
  }

  if ( transformRepository == NULL )
  {
    LOG_ERROR("Failed to add tracked frame to volume - input transform repository is NULL"); 
    return PLUS_FAIL; 
  }

  bool isMatrixValid(false); 
  vtkSmartPointer<vtkMatrix4x4> imageToReferenceTransformMatrix=vtkSmartPointer<vtkMatrix4x4>::New();
  if ( transformRepository->GetTransform(imageToReferenceTransformName, imageToReferenceTransformMatrix, &isMatrixValid ) != PLUS_SUCCESS )
  {
    std::string strImageToReferenceTransformName; 
    imageToReferenceTransformName.GetTransformName(strImageToReferenceTransformName); 
    LOG_ERROR("Failed to get transform '"<<strImageToReferenceTransformName<<"' from transform repository"); 
    return PLUS_FAIL; 
  }

  if ( insertedIntoVolume != NULL )
  {
    *insertedIntoVolume = isMatrixValid; 
  }

  if ( !isMatrixValid )
  {
    // Insert only valid frame into volume
    std::string strImageToReferenceTransformName; 
    imageToReferenceTransformName.GetTransformName(strImageToReferenceTransformName); 
    LOG_DEBUG("Transform '"<<strImageToReferenceTransformName<<"' is invalid for the current frame, therefore this frame is not be inserted into the volume"); 
    return PLUS_SUCCESS; 
  }

  vtkImageData* frameImage=frame->GetImageData()->GetImage();

  this->Modified();

    // initialize variables 
  double currentLeftAngle = -45;
  double currentRightAngle = 45;
  double detectedLeftAngle = 0;
  double detectedRightAngle = 0;
  double xOrigin = 0;
  double yOrigin = 0;
  
  // std::cout << "Initialized some variables" << std::endl;

  // get current fan angle and origin
  this->Reconstructor->GetFanAnglesDeg(currentLeftAngle, currentRightAngle);
  this->Reconstructor->GetFanOrigin(xOrigin, yOrigin);
  
  // std::cout << "Initial Left Angle" << currentLeftAngle << std::endl;
  // std::cout << "Initial Right Angle" << currentRightAngle << std::endl;
  // solicit new fan angles
  this->FanAngleDetector(frameImage, currentLeftAngle, currentRightAngle, detectedLeftAngle, detectedRightAngle, xOrigin, yOrigin);

  std::cout << "Final Left Angle" << detectedLeftAngle << std::endl;
  std::cout << "Final Right Angle" << detectedRightAngle << std::endl;  

  std::cout << "We finished the voodoo?" << std::endl;

  // implement and wrap up 
  this->Reconstructor->SetFanAnglesDeg(detectedLeftAngle, detectedRightAngle);
  PlusStatus temp;
  
  if (detectedLeftAngle > -10.0 && detectedRightAngle < 10.0) {

	temp = PLUS_SUCCESS;
	std::cout << ">>>>>>>>>>>>>BLANK FRAME SKIPPED<<<<<<<<<<<<<<<<<" << std::endl;

  } else {
  
	temp = this->Reconstructor->InsertSlice(frameImage, imageToReferenceTransformMatrix);
  
  }  

  this->Reconstructor->SetFanAnglesDeg(currentLeftAngle, currentRightAngle); 
  return temp;

  // restore original fan angles
}

//----------------------------------------------------------------------------
PlusStatus vtkVolumeReconstructor::UpdateReconstructedVolume()
{
  // Load reconstructed volume if the algorithm configuration was modified since the last load
  //   MTime is updated whenever a new frame is added or configuration modified
  //   ReconstructedVolumeUpdatedTime is updated whenever a reconstruction was completed
  if (this->ReconstructedVolumeUpdatedTime >= this->GetMTime())
  {
    // reconstruction is already up-to-date
    return PLUS_SUCCESS;
  }

  if (this->FillHoles)
  {
    if (this->GenerateHoleFilledVolume() != PLUS_SUCCESS)
    {
      LOG_ERROR("Failed to generate hole filled volume!");
      return PLUS_FAIL;
    }
  }
  else
  {
    this->ReconstructedVolume->DeepCopy(this->Reconstructor->GetReconstructedVolume());
  }

  this->ReconstructedVolumeUpdatedTime = this->GetMTime();

  return PLUS_SUCCESS; 
}

//----------------------------------------------------------------------------
PlusStatus vtkVolumeReconstructor::GetReconstructedVolume(vtkImageData* volume)
{
  if (this->UpdateReconstructedVolume() != PLUS_SUCCESS)
  {
    LOG_ERROR("Failed to load reconstructed volume");
    return PLUS_FAIL;
  }

  volume->DeepCopy(this->ReconstructedVolume);

  return PLUS_SUCCESS; 
}

//----------------------------------------------------------------------------
PlusStatus vtkVolumeReconstructor::GenerateHoleFilledVolume()
{
  LOG_INFO("Hole Filling has begun");
  this->HoleFiller->SetReconstructedVolume(this->Reconstructor->GetReconstructedVolume());
  this->HoleFiller->SetAccumulationBuffer(this->Reconstructor->GetAccumulationBuffer());
  this->HoleFiller->Update();
  LOG_INFO("Hole Filling has finished");

  this->ReconstructedVolume->DeepCopy(HoleFiller->GetOutput());

  return PLUS_SUCCESS; 
}

//----------------------------------------------------------------------------
PlusStatus vtkVolumeReconstructor::ExtractGrayLevels(vtkImageData* reconstructedVolume)
{  
  if (this->UpdateReconstructedVolume() != PLUS_SUCCESS)
  {
    LOG_ERROR("Failed to load reconstructed volume");
    return PLUS_FAIL;
  }

  vtkSmartPointer<vtkImageExtractComponents> extract = vtkSmartPointer<vtkImageExtractComponents>::New();          

  // Keep only first component (the other component is the alpha channel)
  extract->SetComponents(0);
  extract->SetInputData_vtk5compatible(this->ReconstructedVolume);
  extract->Update();

  reconstructedVolume->DeepCopy(extract->GetOutput());

  return PLUS_SUCCESS;
}

//----------------------------------------------------------------------------
PlusStatus vtkVolumeReconstructor::ExtractAlpha(vtkImageData* reconstructedVolume)
{
  if (this->UpdateReconstructedVolume() != PLUS_SUCCESS)
  {
    LOG_ERROR("Failed to load reconstructed volume");
    return PLUS_FAIL;
  }

  vtkSmartPointer<vtkImageExtractComponents> extract = vtkSmartPointer<vtkImageExtractComponents>::New();          

  // Extract the second component (the alpha channel)
  extract->SetComponents(1);
  extract->SetInputData_vtk5compatible(this->ReconstructedVolume);
  extract->Update();

  reconstructedVolume->DeepCopy(extract->GetOutput());

  return PLUS_SUCCESS;
}

//----------------------------------------------------------------------------
PlusStatus vtkVolumeReconstructor::SaveReconstructedVolumeToMetafile(const char* filename, bool alpha/*=false*/, bool useCompression/*=true*/)
{
  vtkSmartPointer<vtkImageData> volumeToSave = vtkSmartPointer<vtkImageData>::New();
  if (alpha)
  {
    if (this->ExtractAlpha(volumeToSave) != PLUS_SUCCESS)
    {
      LOG_ERROR("Extracting alpha channel failed!");
      return PLUS_FAIL;
    }
  }
  else
  {
    if (this->ExtractGrayLevels(volumeToSave) != PLUS_SUCCESS)
    {
      LOG_ERROR("Extracting alpha channel failed!");
      return PLUS_FAIL;
    }
  }
  return SaveReconstructedVolumeToMetafile(volumeToSave, filename, useCompression);
}

//----------------------------------------------------------------------------
PlusStatus vtkVolumeReconstructor::SaveReconstructedVolumeToMetafile(vtkImageData* volumeToSave, const char* filename, bool useCompression/*=true*/)
{
  if (volumeToSave==NULL)
  {
    LOG_ERROR("vtkVolumeReconstructor::SaveReconstructedVolumeToMetafile: invalid input image");
    return PLUS_FAIL;
  }

  MET_ValueEnumType scalarType = MET_NONE;
  switch (volumeToSave->GetScalarType())
  {
  case VTK_UNSIGNED_CHAR: scalarType = MET_UCHAR; break;
  case VTK_FLOAT: scalarType = MET_FLOAT; break;
  default:
    LOG_ERROR("Scalar type is not supported!");
    return PLUS_FAIL;
  }

  MetaImage metaImage(volumeToSave->GetDimensions()[0], volumeToSave->GetDimensions()[1], volumeToSave->GetDimensions()[2],
                                       volumeToSave->GetSpacing()[0], volumeToSave->GetSpacing()[1], volumeToSave->GetSpacing()[2],
                                       scalarType, 1, volumeToSave->GetScalarPointer());
  metaImage.Origin(volumeToSave->GetOrigin());
  // By definition, LPS orientation in DICOM sense = RAI orientation in MetaIO. See details at:
  // http://www.itk.org/Wiki/Proposals:Orientation#Some_notes_on_the_DICOM_convention_and_current_ITK_usage
  metaImage.AnatomicalOrientation("RAI");
  metaImage.BinaryData(true);
  metaImage.CompressedData(useCompression);
  metaImage.ElementDataFileName("LOCAL");
  if (metaImage.Write(filename) == false)
  {
    LOG_ERROR("Failed to save reconstructed volume in sequence metafile!");
    return PLUS_FAIL;
  }

  return PLUS_SUCCESS;
}

//----------------------------------------------------------------------------
PlusStatus vtkVolumeReconstructor::SaveReconstructedVolumeToVtkFile(const char* filename, bool alpha/*=false*/)
{
  vtkSmartPointer<vtkImageData> volumeToSave = vtkSmartPointer<vtkImageData>::New();

  if (alpha)
  {
    if (this->ExtractAlpha(volumeToSave) != PLUS_SUCCESS)
    {
      LOG_ERROR("Extracting alpha channel failed!");
      return PLUS_FAIL;
    }
  }
  else
  {
    if (this->ExtractGrayLevels(volumeToSave) != PLUS_SUCCESS)
    {
      LOG_ERROR("Extracting alpha channel failed!");
      return PLUS_FAIL;
    }
  }

  vtkSmartPointer<vtkDataSetWriter> writer = vtkSmartPointer<vtkDataSetWriter>::New();
  writer->SetFileTypeToBinary();
  writer->SetInputData_vtk5compatible(volumeToSave);
  writer->SetFileName(filename);
  writer->Update();

  return PLUS_SUCCESS;
}

//----------------------------------------------------------------------------
int* vtkVolumeReconstructor::GetClipRectangleOrigin()
{
  return this->Reconstructor->GetClipRectangleOrigin();
}

//----------------------------------------------------------------------------
int* vtkVolumeReconstructor::GetClipRectangleSize()
{
  return this->Reconstructor->GetClipRectangleSize();
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::Reset()
{
  this->Reconstructor->ResetOutput();
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::SetOutputOrigin(double* origin)
{
  this->Reconstructor->SetOutputOrigin(origin);
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::SetOutputSpacing(double* spacing)
{
  this->Reconstructor->SetOutputSpacing(spacing);
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::SetOutputExtent(int* extent)
{
  this->Reconstructor->SetOutputExtent(extent);
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::SetNumberOfThreads(int numberOfThreads)
{
  this->Reconstructor->SetNumberOfThreads(numberOfThreads);
  this->HoleFiller->SetNumberOfThreads(numberOfThreads);
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::SetClipRectangleOrigin(int* origin)
{
  this->Reconstructor->SetClipRectangleOrigin(origin);
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::SetClipRectangleSize(int* size)
{
  this->Reconstructor->SetClipRectangleSize(size);
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::SetFanAnglesDeg(double* anglesDeg)
{
  this->Reconstructor->SetFanAnglesDeg(anglesDeg);
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::SetFanAngles(double* anglesDeg)
{
  LOG_WARNING("FanAngles volume reconstructor parameter is deprecated. Use FanAnglesDeg instead (with the same value).");
  this->Reconstructor->SetFanAnglesDeg(anglesDeg);
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::SetFanOriginPixel(double* originPixel)
{
  // Image coordinate system has unit spacing, so we can set the pixel values directly
  this->Reconstructor->SetFanOrigin(originPixel);
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::SetFanOrigin(double* originPixel)
{
  LOG_WARNING("FanOrigin volume reconstructor parameter is deprecated. Use FanOriginPixels instead (with the same value).");
  SetFanOriginPixel(originPixel);
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::SetFanRadiusStartPixel(double radiusPixel)
{
  // Image coordinate system has unit spacing, so we can set the pixel values directly
  this->Reconstructor->SetFanRadiusStart(radiusPixel);
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::SetFanRadiusStopPixel(double radiusPixel)
{
  // Image coordinate system has unit spacing, so we can set the pixel values directly
  this->Reconstructor->SetFanRadiusStop(radiusPixel);
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::SetFanDepth(double fanDepthPixel)
{
  LOG_WARNING("FanDepth volume reconstructor parameter is deprecated. Use FanRadiusStopPixels instead (with the same value).");
  SetFanRadiusStopPixel(fanDepthPixel);
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::SetInterpolation(vtkPasteSliceIntoVolume::InterpolationType interpolation)
{
  this->Reconstructor->SetInterpolationMode(interpolation);
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::SetCalculation(vtkPasteSliceIntoVolume::CalculationType calculation)
{
  this->Reconstructor->SetCalculationMode(calculation);
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::SetOptimization(vtkPasteSliceIntoVolume::OptimizationType optimization)
{
  this->Reconstructor->SetOptimization(optimization);
}

//----------------------------------------------------------------------------
void vtkVolumeReconstructor::SetCompounding(bool enable)
{
  this->Reconstructor->SetCompounding(enable?1:0);
}

void vtkVolumeReconstructor::FanAngleDetector(vtkImageData* frameImage, double &presetAngleLeft, double &presetAngleRight, double &outputAngleLeft, double &outputAngleRight, double &xOrigin, double &yOrigin)
{

	int nRadii = 3;
	double testRadius1 = 320;
	double testRadius2 = 380;
	double testRadius3 = 440; //480
	double angularResolution = 1;
	double thetaIncrement1 = angularResolution/testRadius1 * (180 / vtkMath::Pi()); 
	double thetaIncrement2 = angularResolution/testRadius2 * (180 / vtkMath::Pi()); 
	double thetaIncrement3 = angularResolution/testRadius3 * (180 / vtkMath::Pi()); 

	// Angle Detection Parameters
	double nPix  = 20;
	double threshold = 33; //26
	double buffer = 0; //degrees

	// Indecies and definitions 
	int left = 0;
	int right = 1;
	int test1 = 0;
	int test2 = 1;
	int test3 = 2;
	int bad = 0;
	int good = 1;

	// create and initialize vectors of test thetas
	std::vector<double> testTheta1;
	std::vector<double> testTheta2;
	std::vector<double> testTheta3;
	testTheta1.push_back(presetAngleLeft);
	testTheta2.push_back(presetAngleLeft);
	testTheta3.push_back(presetAngleLeft);

	// initialize test (x,y) coordinates based on the first test angle
	int testX1 = vtkMath::Round(xOrigin + testRadius1*sin(vtkMath::RadiansFromDegrees(testTheta1[0])));
	int testX2 = vtkMath::Round(xOrigin + testRadius2*sin(vtkMath::RadiansFromDegrees(testTheta2[0])));
	int testX3 = vtkMath::Round(xOrigin + testRadius3*sin(vtkMath::RadiansFromDegrees(testTheta3[0])));
	int testY1 = vtkMath::Round(yOrigin + testRadius1*cos(vtkMath::RadiansFromDegrees(testTheta1[0])));
	int testY2 = vtkMath::Round(yOrigin + testRadius2*cos(vtkMath::RadiansFromDegrees(testTheta2[0])));
	int testY3 = vtkMath::Round(yOrigin + testRadius3*cos(vtkMath::RadiansFromDegrees(testTheta3[0])));

	// print screens for debugging
	// std::cout << testX1 << std::endl;
	// std::cout << testY1 << std::endl;
	// std::cout << frameImage->GetScalarTypeAsString() << std::endl;

	// create and initialize testValue vectors
	std::vector<unsigned char> testValue1; 
	std::vector<unsigned char> testValue2;
	std::vector<unsigned char> testValue3;
	testValue1.push_back(*static_cast<unsigned char*>(frameImage->GetScalarPointer(testX1,testY1,0)));
	testValue2.push_back(*static_cast<unsigned char*>(frameImage->GetScalarPointer(testX2,testY2,0)));
	testValue3.push_back(*static_cast<unsigned char*>(frameImage->GetScalarPointer(testX2,testY2,0)));

	// run three while loops to populate the testValue vectors. Note that test x's and y's are not being retained.
	int i = 0;

	while (testTheta1[i] + thetaIncrement1 <= presetAngleRight) {

		testTheta1.push_back(testTheta1[i] + thetaIncrement1);
		testX1 = vtkMath::Round(xOrigin + testRadius1*sin(vtkMath::RadiansFromDegrees(testTheta1[i+1])));
		testY1 = vtkMath::Round(yOrigin + testRadius1*cos(vtkMath::RadiansFromDegrees(testTheta1[i+1])));
		testValue1.push_back(*static_cast<unsigned char*>(frameImage->GetScalarPointer(testX1,testY1,0)));
		//std::cout << (int) testValue1[i] << std::endl;
		i++;

	}

	i = 0;

	while (testTheta2[i] + thetaIncrement2 <= presetAngleRight) {

		testTheta2.push_back(testTheta2[i] + thetaIncrement2);
		testX2 = vtkMath::Round(xOrigin + testRadius2*sin(vtkMath::RadiansFromDegrees(testTheta2[i+1])));
		testY2 = vtkMath::Round(yOrigin + testRadius2*cos(vtkMath::RadiansFromDegrees(testTheta2[i+1])));
		testValue2.push_back(*static_cast<unsigned char*>(frameImage->GetScalarPointer(testX2,testY2,0)));
		i++;

	}

	i = 0;

	while (testTheta3[i] + thetaIncrement3 <= presetAngleRight) {

		testTheta3.push_back(testTheta3[i] + thetaIncrement3);
		testX3 = vtkMath::Round(xOrigin + testRadius3*sin(vtkMath::RadiansFromDegrees(testTheta3[i+1])));
		testY3 = vtkMath::Round(yOrigin + testRadius3*cos(vtkMath::RadiansFromDegrees(testTheta3[i+1])));
		testValue3.push_back(*static_cast<unsigned char*>(frameImage->GetScalarPointer(testX3,testY3,0)));
		i++;

	}

	// Print lines for debugging
	// std::cout << testTheta1.size() << std::endl;
	// std::cout << testTheta2.size() << std::endl;
	// std::cout << testTheta3.size() << std::endl;
	// Find the limiting fan angles

	// create fan angles vector and initialize values to zero
	std::vector< std::vector<double> > fanAngles(nRadii, std::vector<double>(2));

	for (int i = 0; i < nRadii; i++) {

		for (int j = 0; j < 2; j++) {

			fanAngles[i][j] = 0;

		}
	}

	// variables for use in the loop
	int nTheta1 = testTheta1.size();
	int nTheta2 = testTheta2.size();
	int nTheta3 = testTheta3.size();

	// create and initialize state variable
	int state1 = bad;
	int state2 = bad;
	int state3 = bad;

	// create and initialize angle log
	std::vector<double> left_log_1;
	std::vector<double> right_log_1;
	std::vector<double> left_log_2;
	std::vector<double> right_log_2;
	std::vector<double> left_log_3;
	std::vector<double> right_log_3;

	left_log_1.push_back(presetAngleRight);
	left_log_2.push_back(presetAngleRight);
	left_log_3.push_back(presetAngleRight);
	right_log_1.push_back(presetAngleRight);
	right_log_2.push_back(presetAngleRight);
	right_log_3.push_back(presetAngleRight);

	int j = nPix;

	while (true) {

		if (j <= nTheta1) {

			double testCount = 0;

			for (int k = 0; k < nPix; k++) {

				if (testValue1[j - nPix + k] >= threshold) {testCount++;}

			}

			if (state1 == bad && testCount >= nPix/2.0) {

				left_log_1[left_log_1.size() - 1] = testTheta1[j - nPix];
				state1 = good;

			} else if (state1 == good && testCount <= nPix/2.0 ) {

				right_log_1[right_log_1.size() - 1] = testTheta1[j - 1];
				left_log_1.push_back(presetAngleRight);
				right_log_1.push_back(presetAngleRight);
				state1 = bad;

			}



		}

		if (j <= nTheta2) {

			double testCount = 0;

			for (int k = 0; k < nPix; k++) {

				if (testValue2[j - nPix + k] >= threshold) {testCount++;}

			}

			if (state2 == bad && testCount >= nPix/2.0) {

				left_log_2[left_log_2.size() - 1] = testTheta2[j - nPix];
				state2 = good;

			} else if (state2 == good && testCount <= nPix/2.0 ) {

				right_log_2[right_log_2.size() - 1] = testTheta2[j - 1];
				left_log_2.push_back(presetAngleRight);
				right_log_2.push_back(presetAngleRight);
				state2 = bad;

			}



		}

		if (j <= nTheta3) {

			double testCount = 0;

			for (int k = 0; k < nPix; k++) {

				if (testValue3[j - nPix + k] >= threshold) {testCount++;}

			}

			if (state3 == bad && testCount >= nPix/2.0) {

				left_log_3[left_log_3.size() - 1] = testTheta3[j - nPix];
				state3 = good;

			} else if (state3 == good && testCount <= nPix/2.0 ) {

				right_log_3[right_log_3.size() - 1] = testTheta3[j - 1];
				left_log_3.push_back(presetAngleRight);
				right_log_3.push_back(presetAngleRight);
				state3 = bad;

			}


		} else {

			break;

		}

		j++; 

	}

	std::cout << "Completed the detection loop?!" << std::endl;

	// pick the largest continuous area of pixel value for each test radius and send it to fanAngles

	fanAngles[test1][left] = left_log_1[0];
	fanAngles[test1][right] = right_log_1[0];

	for (int k = 0; k < left_log_1.size(); k++) {

		if ( (right_log_1[k] - left_log_1[k]) > (fanAngles[test1][right] - fanAngles[test1][left]) ) {

			fanAngles[test1][left] = left_log_1[k];
			fanAngles[test1][right] = right_log_1[k];

		}

	}

	// std::cout << "log 1" << std::endl;
	// std::cout << left_log_1 << std::endl;
	// std::cout << right_log_1 << std::endl;

	fanAngles[test2][left] = left_log_2[0];
	fanAngles[test2][right] = right_log_2[0];

	for (int k = 0; k < left_log_2.size(); k++) {

		if ( (right_log_2[k] - left_log_2[k]) > (fanAngles[test2][right] - fanAngles[test2][left]) ) {

			fanAngles[test2][left] = left_log_2[k];
			fanAngles[test2][right] = right_log_2[k];

		}

	}

	// std::cout << "log 2" << std::endl;
	// std::cout << left_log_2 << std::endl;
	// std::cout << right_log_2 << std::endl;

	fanAngles[test3][left] = left_log_3[0];
	fanAngles[test3][right] = right_log_3[0];

	for (int k = 0; k < left_log_3.size(); k++) {

		if ( (right_log_3[k] - left_log_3[k]) > (fanAngles[test3][right] - fanAngles[test3][left]) ) {

			fanAngles[test3][left] = left_log_3[k];
			fanAngles[test3][right] = right_log_3[k];

		}

	}


	// std::cout << "log 3" << std::endl;
	// std::cout << left_log_3 << std::endl;
	// std::cout << right_log_3 << std::endl;

	// search for the broadest of the three angular windows and set them as the outputs
	// broadest in case there's a ventricle or something
	outputAngleLeft = fanAngles[test1][left];
	outputAngleRight = fanAngles[test1][right];

	for (int ii = 1; ii < nRadii; ii++) {

		if ( (fanAngles[ii][right] - fanAngles[ii][left]) < (outputAngleRight - outputAngleLeft) ) {

			outputAngleLeft = fanAngles[ii][left];
			outputAngleRight = fanAngles[ii][right];

		}

	}

	// set the blank frame case (i.e. both angles = presetRightAngle) to be 0 degrees for consistency
	if (outputAngleLeft == presetAngleRight && outputAngleRight == presetAngleRight) {

		outputAngleLeft = 0;
		outputAngleRight = 0;

	}

	// apply a buffer if necessary 
	if (outputAngleLeft != 0) {outputAngleLeft = outputAngleLeft - buffer;}
	if (outputAngleRight != 0) {outputAngleRight = outputAngleRight + buffer;}

}
