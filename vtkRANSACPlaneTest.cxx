#include <vtkSmartPointer.h>

#include "vtkRANSACPlane.h"

int main (int argc, char *argv[])
{
  
  vtkSmartPointer<vtkRANSACPlane> RANSACPlane = vtkSmartPointer<vtkRANSACPlane>::New();
  
  if(RANSACPlane->GetInlierThreshold() != 1.0)
    {
    return EXIT_FAILURE;
    }
  
  if(RANSACPlane->GetMaxIterations() != 1000)
    {
    return EXIT_FAILURE;
    }

  if(RANSACPlane->GetGoodEnough() != 1.0)
    {
    return EXIT_FAILURE;
    }
    
  return EXIT_SUCCESS;
  
  return 0;
}
