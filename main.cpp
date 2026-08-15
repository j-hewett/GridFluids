#include "GridWidget.h"
#include "Fluid.h"
#include <QApplication>
#include <QTimer>


int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    int gridSize = 40;

    GridWidget grid(gridSize, gridSize);
    grid.setWindowTitle("Fluid");
    grid.resize(800,800);
    grid.show();

    Fluid fluid(gridSize, 5000);

    QTimer timer;
    float dt = 1.0f/60.0f;

    QObject::connect(&timer, &QTimer::timeout, [&]()
    {
        auto cellType = fluid.findFluid(dt);
        std::vector<int> occupancy(cellType.size());
        for (size_t i = 0; i < cellType.size(); ++i)
            occupancy[i] = (cellType[i] == Fluid::FLUID_CELL) ? 1 : 0;
        grid.setCells(occupancy);
    });

    timer.start(1000*dt);

    return app.exec();
}
