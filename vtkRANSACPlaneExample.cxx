#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

#include "vtkRANSACPlane.h"

int main (int argc, char *argv[])
{
  //verify command line arguments
  if(argc != 3)
    {
    cout << "Required arguments: InputFilename OutputFilename" << endl;
    return EXIT_FAILURE;
    }
  
  //parse command line arguments
  vtkstd::string inputFilename = argv[1];
  vtkstd::string outputFilename = argv[2];
  
  //read the input file
  vtkSmartPointer<vtkXMLPolyDataReader> reader = 
      vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(inputFilename.c_str());
  reader->Update();
  
  //vtkPolyData* inputPoints = Reader->GetOutput();
  
  //estimate normals
  vtkSmartPointer<vtkRANSACPlane> rANSACPlane = 
      vtkSmartPointer<vtkRANSACPlane>::New();
  rANSACPlane->SetInlierThreshold(0.1);
  rANSACPlane->SetMaxIterations(1000);
  //rANSACPlane->SetInput(inputPoints);
  rANSACPlane->SetInputConnection(reader->GetOutputPort());
  rANSACPlane->Update();
  
  
  //write the output file with the estimated normals
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = 
      vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(outputFilename.c_str());
  writer->SetInput(rANSACPlane->GetOutput());
  writer->Write();
  
  return EXIT_SUCCESS;
}
