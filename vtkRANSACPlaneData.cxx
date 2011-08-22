#include <vtkPolyData.h>
#include <vtkPlaneSource.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkMath.h>
#include <vtkCellArray.h>

int main (int argc, char *argv[])
{
  //verify command line arguments
  if(argc != 2)
    {
    cout << "Required arguments: OutputFilename" << endl;
    return EXIT_FAILURE;
    }
  
  //parse command line arguments
  vtkstd::string outputFilename = argv[1];
  
  vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
  planeSource->SetCenter(1.0, 0.0, 0.0);
  planeSource->SetNormal(1.0, 0.0, 0.0);
  unsigned int res = 10;
  planeSource->SetResolution(res, res);
  planeSource->Update();
  
  vtkPolyData* oldPolydata = planeSource->GetOutput();
  vtkPoints* oldPoints = oldPolydata->GetPoints();
  vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
  
  //add erratic points and noise
  double noise = .1;
  for(unsigned int i = 0; i < oldPoints->GetNumberOfPoints(); i++)
    {
    double p[3];
    oldPoints->GetPoint(i, p);
    for(unsigned int j = 0; j < 3; j++)
      {
      p[j] = p[j] + vtkMath::Random(-noise, noise);
      }
      
      vtkIdType pid[1];
      pid[0] = newPoints->InsertNextPoint(p[0], p[1], p[2]);
      vertices->InsertNextCell ( 1,pid );
    }
    
  unsigned int numOutliers = oldPoints->GetNumberOfPoints() * .1;
  for(unsigned int i = 0; i < numOutliers; i++)
    {
    double p[3];
    for(unsigned int j = 0; j < 3; j++)
      {
      p[j] = vtkMath::Random(-5, 5);
      }
      vtkIdType pid[1];
      pid[0] = newPoints->InsertNextPoint(p[0], p[1], p[2]);
      vertices->InsertNextCell ( 1,pid );
    }
  
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(newPoints);
  polydata->SetVerts(vertices);
    
  //write the output file with the estimated normals
  vtkSmartPointer<vtkXMLPolyDataWriter> Writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  Writer->SetFileName(outputFilename.c_str());
  Writer->SetInput(polydata);
  Writer->Write();
  
  return 0;
}
