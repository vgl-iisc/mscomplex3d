#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkRendererCollection.h>
#include <vtkDataSetMapper.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIdTypeArray.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCommand.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPlaneSource.h>
#include <vtkCellPicker.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkProperty.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkExtractSelection.h>
#include <vtkObjectFactory.h>
#include <vtkInteractorStyleRubberBandPick.h>
#include <vtkRenderedAreaPicker.h>
#include "vtkCallbackCommand.h"
#include <vtkHardwareSelector.h>
#include "vtkExtractSelectedPolyDataIds.h"
#include "vtkSetGet.h"

#define MY_CREATE_NEW(class, variable)\
  vtkSmartPointer<class> variable = vtkSmartPointer<class>::New();

#define MY_NEW(class)\
  vtkSmartPointer<class>::New();


struct pers_diagram
{
  vtkSmartPointer<vtkPoints>          points;
  vtkSmartPointer<vtkPolyData>        polydata;
  vtkSmartPointer<vtkPolyDataMapper>  mapper;
  vtkSmartPointer<vtkActor>           actor;
  vtkSmartPointer<vtkRenderer>        renderer;

//  vtkSmartPointer<vtkPoints>          points;
//  vtkSmartPointer<vtkPolyData>        polydata;
  vtkSmartPointer<vtkPolyDataMapper>  sel_mapper;
  vtkSmartPointer<vtkActor>           sel_actor;


  void build_pipeline()
  {
   // Create the geometry of a point (the coordinate)
    points   = MY_NEW(vtkPoints);
    polydata = MY_NEW(vtkPolyData);
    mapper   = MY_NEW(vtkPolyDataMapper);
    actor    = MY_NEW(vtkActor);
    renderer = MY_NEW(vtkRenderer);

    const float p[] = {
      0,0,1,
      0,1,1,
      1,0,1,
      1,1,1,
    };


    // Create the topology of the point (a vertex)

    for(int i = 0 ; i < sizeof(p)/(3*sizeof(float));++i)
    {
      points->InsertNextPoint(&p[3*i]);
    }

    // Create a polydata object
    polydata->SetPoints(points);

    vtkIdType quad[] = {0,1,2,3};

    polydata->Allocate(100);

    polydata->InsertNextCell(VTK_VERTEX,1,quad);
    polydata->InsertNextCell(VTK_VERTEX,1,quad+1);
    polydata->InsertNextCell(VTK_VERTEX,1,quad+2);
    polydata->InsertNextCell(VTK_VERTEX,1,quad+3);

    polydata->InsertNextCell(VTK_TRIANGLE,3,quad);
    polydata->InsertNextCell(VTK_TRIANGLE,3,quad+1);

    // Visualize
#if VTK_MAJOR_VERSION <= 5
    mapper->SetInput(polydata);
#else
    mapper->SetInputData(point);
#endif

    actor->SetMapper(mapper);
    actor->GetProperty()->SetPointSize(8);
    actor->GetProperty()->SetEdgeVisibility(true);

    renderer->AddActor(actor);
    renderer->ResetCamera();
    renderer->SetBackground(0,0,0); // Blue
  }

  void EndPick()
  {

    MY_CREATE_NEW(vtkHardwareSelector, sel);
    sel->SetRenderer(renderer);

    double x0 = renderer->GetPickX1();
    double y0 = renderer->GetPickY1();
    double x1 = renderer->GetPickX2();
    double y1 = renderer->GetPickY2();

    sel->SetArea(static_cast<int>(x0),
                 static_cast<int>(y0),
                 static_cast<int>(x1),
                 static_cast<int>(y1));

    vtkSmartPointer<vtkSelection> res;
    res.TakeReference(sel->Select());
    if (!res)
    {
      cerr << "Selection not supported." << endl;
      return;
    }

    cerr << "x0 " << x0 << " y0 " << y0 << "\t";
    cerr << "x1 " << x1 << " y1 " << y1 << endl;

    vtkSelectionNode *cellids = res->GetNode(0);
    MY_CREATE_NEW(vtkExtractSelectedPolyDataIds, extr);

    if (cellids)
    {
      vtkIdTypeArray *id_arr =
          dynamic_cast<vtkIdTypeArray*>(cellids->GetSelectionList());

      for( int i = 0 ; i < id_arr->GetDataSize(); ++i)
      {
        vtkIdType id = id_arr->GetValue(i);
        int     type = polydata->GetCellType(id);
        cout<<"id:type = "<<id<<":" << type <<endl;
      }
    }
  }

  inline vtkSmartPointer<vtkActor>    get_actor(){return actor;}
  inline vtkSmartPointer<vtkPolyData> get_polydata(){return polydata;}
  inline vtkSmartPointer<vtkRenderer> get_renderer(){return renderer;}
};

void EndPick (vtkObject *vtkNotUsed( caller ),
                     unsigned long vtkNotUsed(eventId),
                     void * clientData, void *)
{
  pers_diagram *ptr = (pers_diagram*) clientData;
  ptr->EndPick();
}

int main_old (int, char *[])
{

  pers_diagram pd;
  pd.build_pipeline();

  MY_CREATE_NEW(vtkRenderWindow,renderWindow);
  MY_CREATE_NEW(vtkRenderWindowInteractor,renderWindowInteractor);
  MY_CREATE_NEW(vtkInteractorStyleRubberBandPick,rbp);
  MY_CREATE_NEW(vtkRenderedAreaPicker,areaPicker);

  renderWindow->AddRenderer(pd.get_renderer());
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->Initialize();

  // Set the custom stype to use for interaction.
  renderWindowInteractor->SetInteractorStyle(rbp);
  renderWindowInteractor->SetPicker(areaPicker);


  //pass pick events to the HardwareSelector
  MY_CREATE_NEW(vtkCallbackCommand, cbc);
  cbc->SetCallback(EndPick);
  cbc->SetClientData(&pd);
  renderWindowInteractor->AddObserver(vtkCommand::EndPickEvent,cbc);

  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}

#include <QApplication>
#include "mainwindow.h"

int main( int argc, char** argv )
{
  // QT Stuff
  QApplication app( argc, argv );

  MainWindow mw;
  mw.show();

  return app.exec();
}

