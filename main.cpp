#include "GridWidget.h"
#include "Fluid.h"
#include <QApplication>
#include <QTimer>


int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    int gridSize = 30;

    GridWidget grid(gridSize, gridSize);
    grid.setWindowTitle("Fluid");
    grid.resize(800,800);
    grid.show();

    Fluid fluid(gridSize, 2);

    QTimer timer;
    float dt = 0.1;

    QObject::connect(&timer, &QTimer::timeout, [&]()
    {
        auto gridState = fluid.findFluid(dt);
        grid.setCells(gridState);
    });

    timer.start(1000*dt);

    return app.exec();
}
