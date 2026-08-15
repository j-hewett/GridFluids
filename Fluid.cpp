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

    //u component: staggered offset (0, 0.5*cellSize)
    for (const auto& p : particles)
    {
        float gx = p.pos.x / cellSize;              // no offset in x
        float gy = p.pos.y / cellSize - 0.5f;        // offset by half a cell in y

        int i0 = static_cast<int>(std::floor(gx));
        int j0 = static_cast<int>(std::floor(gy));
        int i1 = i0 + 1;
        int j1 = j0 + 1;

        float tx = gx - i0;
        float ty = gy - j0;

        // clamp to valid u-grid index range: i in [0,size], j in [0,size-1]
        i0 = clampInt(i0, 0, size);
        i1 = clampInt(i1, 0, size);
        j0 = clampInt(j0, 0, size - 1);
        j1 = clampInt(j1, 0, size - 1);

        float w00 = (1 - tx) * (1 - ty);
        float w10 = tx * (1 - ty);
        float w01 = (1 - tx) * ty;
        float w11 = tx * ty;

        velocitiesX[idxX(i0, j0)] += p.vel.x * w00;  weightsX[idxX(i0, j0)] += w00;
        velocitiesX[idxX(i1, j0)] += p.vel.x * w10;  weightsX[idxX(i1, j0)] += w10;
        velocitiesX[idxX(i0, j1)] += p.vel.x * w01;  weightsX[idxX(i0, j1)] += w01;
        velocitiesX[idxX(i1, j1)] += p.vel.x * w11;  weightsX[idxX(i1, j1)] += w11;
    }

    //v component: staggered offset (0.5*cellSize, 0)
    for (const auto& p : particles)
    {
        float gx = p.pos.x / cellSize - 0.5f;
        float gy = p.pos.y / cellSize;

        int i0 = static_cast<int>(std::floor(gx));
        int j0 = static_cast<int>(std::floor(gy));
        int i1 = i0 + 1;
        int j1 = j0 + 1;

        float tx = gx - i0;
        float ty = gy - j0;

        //clamp to valid v-grid index range: i in [0,size-1], j in [0,size]
        i0 = clampInt(i0, 0, size - 1);
        i1 = clampInt(i1, 0, size - 1);
        j0 = clampInt(j0, 0, size);
        j1 = clampInt(j1, 0, size);

        float w00 = (1 - tx) * (1 - ty);
        float w10 = tx * (1 - ty);
        float w01 = (1 - tx) * ty;
        float w11 = tx * ty;

        velocitiesY[idxY(i0, j0)] += p.vel.y * w00;  weightsY[idxY(i0, j0)] += w00;
        velocitiesY[idxY(i1, j0)] += p.vel.y * w10;  weightsY[idxY(i1, j0)] += w10;
        velocitiesY[idxY(i0, j1)] += p.vel.y * w01;  weightsY[idxY(i0, j1)] += w01;
        velocitiesY[idxY(i1, j1)] += p.vel.y * w11;  weightsY[idxY(i1, j1)] += w11;
    }

    //normalise by accumulated weight
    for (size_t idx = 0; idx < velocitiesX.size(); ++idx)
        if (weightsX[idx] > 0.0f)
            velocitiesX[idx] /= weightsX[idx];

    for (size_t idx = 0; idx < velocitiesY.size(); ++idx)
        if (weightsY[idx] > 0.0f)
            velocitiesY[idx] /= weightsY[idx];
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
