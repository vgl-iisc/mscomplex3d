#include <vtkDataObjectToTable.h>
#include <vtkElevationFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkQtTableView.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkSphereSource.h>
#include <vtkCubeSource.h>
#include <vtkInteractorStyleRubberBandPick.h>
#include <vtkRenderedAreaPicker.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkPointData.h>

#include <QFileDialog>
#include <QDir>

#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

#include <ui_mainwindow.h>
#include <mainwindow.h>

#define MY_CREATE_NEW(class, variable)\
  vtkSmartPointer<class> variable = vtkSmartPointer<class>::New();

inline std::string msgraphname_to_minmfoldname(std::string fname)
{
  boost::algorithm::replace_last(fname,"graph.bin","min.raw");
  return fname;
}

inline std::string msgraphname_to_maxmfoldname(std::string fname)
{
  boost::algorithm::replace_last(fname,"graph.bin","max.raw");
  return fname;
}

// Constructor
MainWindow::MainWindow()
{
  this->ui = new Ui_MainWindow;
  this->ui->setupUi(this);

  // Set up action signals and slots
  connect(this->ui->actionExit, SIGNAL(triggered()), this, SLOT(slotExit()));

  // cube
  MY_CREATE_NEW(vtkCubeSource,cubeSource);
  MY_CREATE_NEW(vtkPolyDataMapper,cubeMapper);
  MY_CREATE_NEW(vtkActor,cubeActor);
  MY_CREATE_NEW(vtkRenderer,rightRenderer);

  cubeSource->Update();
  cubeMapper->SetInputConnection(cubeSource->GetOutputPort());
  cubeActor->SetMapper(cubeMapper);
  rightRenderer->AddActor(cubeActor);

  this->ui->qvtkWidgetRight->GetRenderWindow()->AddRenderer(rightRenderer);

  // add pick interactor to left renderer
  MY_CREATE_NEW(QVTKInteractor,renderWindowInteractor);
  MY_CREATE_NEW(vtkInteractorStyleRubberBandPick,rbp);
  MY_CREATE_NEW(vtkRenderedAreaPicker,areaPicker);

  renderWindowInteractor->SetRenderWindow
      (this->ui->qvtkWidgetLeft->GetRenderWindow());
  renderWindowInteractor->Initialize();

  // Set the custom stype to use for interaction.
  renderWindowInteractor->SetInteractorStyle(rbp);
  renderWindowInteractor->SetPicker(areaPicker);
};

void MainWindow::slotExit()
{
  qApp->exit();
}

void MainWindow::on_actionOpenMsComplex_triggered(bool )
{
  std::string msc_file = QFileDialog::getOpenFileName
      (this,tr("Select mscomplex file"),QDir::currentPath(),
       "Mscomplex (*.graph.bin)").toStdString();

  if(msc_file == "")
    return;

  using namespace grid;

  m_msc.reset(new grid::mscomplex_t);

  m_msc->load(msc_file);

  MY_CREATE_NEW(vtkPoints,pers_pts);
  MY_CREATE_NEW(vtkUnsignedCharArray,pers_pt_colors);
  MY_CREATE_NEW(vtkPolyData,pers_polydata);
  MY_CREATE_NEW(vtkPolyDataMapper,pers_mapper);
  MY_CREATE_NEW(vtkActor,pers_actor);
  MY_CREATE_NEW(vtkRenderer,pers_renderer);

  pers_pt_colors->SetNumberOfComponents(3);
  pers_pt_colors->SetName("Colors");

  mscomplex_t &msc = *m_msc;
  double pt[] = {0,0,0};

  unsigned char ub_red   [] = {255,0,0};
  unsigned char ub_green [] = {0,255,0};
  unsigned char ub_blue  [] = {0,0,255};

  BOOST_FOREACH(int_pair_t pr, msc.m_canc_list)
  {
    if(msc.index(pr[0]) > msc.index(pr[1]))
      std::swap(pr[0],pr[1]);

    pt[0] = msc.fn(pr[0]);
    pt[1] = msc.fn(pr[1]);

    switch(msc.index(pr[0]))
    {
    case 0:pers_pt_colors->InsertNextTupleValue(ub_red);break;
    case 1:pers_pt_colors->InsertNextTupleValue(ub_green);break;
    case 2:pers_pt_colors->InsertNextTupleValue(ub_blue);break;
    };

    pers_pts->InsertNextPoint(pt);
  }

  pt[0] = msc.fn_min(); pt[1] = pt[0];
  pers_pts->InsertNextPoint(pt);
  pt[1] = msc.fn_max();
  pers_pts->InsertNextPoint(pt);
  pt[0] = pt[1];
  pers_pts->InsertNextPoint(pt);

  pers_polydata->SetPoints(pers_pts);
  pers_polydata->GetPointData()->SetScalars(pers_pt_colors);

  pers_polydata->Allocate(msc.m_canc_list.size() + 4);

  for(vtkIdType i = 0 ; i < msc.m_canc_list.size() ; ++i)
    pers_polydata->InsertNextCell(VTK_VERTEX,1,&i);

  vtkIdType tri[] = {msc.m_canc_list.size(),msc.m_canc_list.size()+1,
                    msc.m_canc_list.size()+2};

  pers_polydata->InsertNextCell(VTK_TRIANGLE,3,tri);

#if VTK_MAJOR_VERSION <= 5
  pers_mapper->SetInput(pers_polydata);
#else
  pers_mapper->SetInputData(pers_pdata);
#endif

  pers_actor->SetMapper(pers_mapper);
  pers_actor->GetProperty()->SetPointSize(3);
  pers_actor->GetProperty()->SetEdgeVisibility(true);
  pers_actor->GetProperty()->SetEdgeColor(1,0,0);
  pers_actor->GetProperty()->SetColor(0,0,1);
  pers_actor->GetProperty()->SetRepresentationToWireframe();

  pers_renderer->AddActor(pers_actor);
  pers_renderer->ResetCamera();
  pers_renderer->SetBackground(1,1,1);

  this->ui->qvtkWidgetLeft->GetRenderWindow()->AddRenderer(pers_renderer);
}

void MainWindow::on_actionOpenVolume_triggered(bool )
{
  if(!m_msc) return;
}
