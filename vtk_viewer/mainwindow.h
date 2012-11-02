#ifndef MainWindow_H
#define MainWindow_H

#include "vtkSmartPointer.h"
#include <QMainWindow>

#include <grid_mscomplex.h>

// Forward Qt class declarations
class Ui_MainWindow;

class MainWindow : public QMainWindow
{
  Q_OBJECT
public:

  // Constructor/Destructor
  MainWindow();
  ~MainWindow() {};

public slots:

  virtual void slotExit();

  void on_actionOpenMsComplex_triggered(bool );

  void on_actionOpenVolume_triggered(bool );

protected:

protected slots:

private:

  // Designer form
  Ui_MainWindow *ui;

  grid::mscomplex_ptr_t m_msc;
};

#endif // MainWindow_H
