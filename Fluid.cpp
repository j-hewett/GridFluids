#include "Fluid.h"
#include <random>

static inline int clampInt(int v, int lo, int hi)
{
    return (v < lo) ? lo : (v > hi) ? hi : v;
}

Fluid::Fluid(int size, int n_particles)
    : size(size),
    cellSize(1.0 / size),
    gridState(size * size, 0.0),
    particles(n_particles),
    velocitiesX((size+1) * size, 0.0),
    velocitiesY(size * (size+1), 0.0),
    pressures(size*size, 0.0),
    weightsX((size+1)*size, 0.0f),
    weightsY(size*(size+1), 0.0f)
{
    genParticles();
}

float Fluid::velDivAtCell(int cellX, int cellY)
{
    float velTop = velocitiesY[idxY(cellX+0, cellY+1)];
    float velLeft = velocitiesX[idxX(cellX+0, cellY+0)];
    float velRight = velocitiesX[idxX(cellX+1, cellY+0)];
    float velBottom = velocitiesY[idxY(cellX+0, cellY+0)];

    float gradX = (velRight - velLeft) / cellSize;
    float gradY = (velTop - velBottom) / cellSize;

    float div = gradX + gradY;
    return div;
}

std::tuple<size_t, size_t> Fluid::findCell(Particle p)
{
    //delicate - explicitly define as tuple
    return {static_cast<size_t>(p.pos.x/cellSize), static_cast<size_t>(p.pos.y/cellSize)};
}

void Fluid::genParticles()
{
    static std::default_random_engine rng(42);
    std::uniform_real_distribution<float> dist(0.0, 1.0);

    for (auto& p : particles)
    {
        p.pos = vec2(dist(rng), dist(rng));
        vec2 randVel = vec2((dist(rng)-0.5)*0.5, (dist(rng)-0.5)*0.5);
        p.vel = randVel;
    }
}

void Fluid::clearGrid()
{
    std::fill(velocitiesX.begin(), velocitiesX.end(), 0.0f);
    std::fill(velocitiesY.begin(), velocitiesY.end(), 0.0f);
    std::fill(weightsX.begin(), weightsX.end(), 0.0f);
    std::fill(weightsY.begin(), weightsY.end(), 0.0f);
}

void Fluid::particles2Grid()
{
    clearGrid();

    for (const auto&p : particles)
    {
        // x
        float gx = p.pos.x / cellSize;
        float gy = (p.pos.y / cellSize) - 0.5f;

        int i = static_cast<int>(gx);
        int j = static_cast<int>(gy);

        //offsets
        float ox = p.pos.x - gx*cellSize;
        float oy = p.pos.y - gy*cellSize;

        float w1 = (1 - (ox/cellSize))*(1 - (oy/cellSize));
        float w2 = ox/cellSize * (1 - (oy/cellSize));
        float w3 = (1 - (ox/cellSize)) * oy/cellSize;
        float w4 = ox/cellSize * oy/cellSize;
    }
}

std::vector<int> Fluid::findFluid(float dt)
{
    std::fill(gridState.begin(), gridState.end(), 0);
    for(auto& p : particles)
    {
        vec2 gravity{0.0f, 1.5f};
        p.integrate(dt, gravity);
        auto coords = findCell(p);
        size_t col = std::get<0>(coords);
        size_t row = std::get<1>(coords);

        if (col < static_cast<size_t>(size) && row < static_cast<size_t>(size))
            gridState[row * size + col] = 1;
    }
    return gridState;
}
