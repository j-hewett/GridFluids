#include "Fluid.h"
#include <random>

static inline int clampInt(int v, int lo, int hi)
{
    return (v < lo) ? lo : (v > hi) ? hi : v;
}

Fluid::Fluid(int size, int n_particles)
    : size(size),
    cellSize(1.0 / size),
    particles(n_particles),
    s(size * size, 1.0f),
    cellType(size * size, AIR_CELL),
    velocitiesX((size+1) * size, 0.0),
    velocitiesY(size * (size+1), 0.0),
    preVelocitiesX((size+1) * size, 0.0f),
    preVelocitiesY(size * (size+1), 0.0f),
    pressures(size*size, 0.0),
    weightsX((size+1)*size, 0.0f),
    weightsY(size*(size+1), 0.0f)
{
    genParticles();
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

void Fluid::initBoundaries()
{
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            bool isWall = (i == 0 || i == size - 1 || j == 0 || j == size - 1);
            s[idxC(i, j)] = isWall ? 0.0f : 1.0f;
        }
    }
}

float Fluid::velDivAtCell(int cellX, int cellY)
{
    if (cellType[idxC(cellX, cellY)] != FLUID_CELL)
        return 0.0f;

    float velTop = velocitiesY[idxY(cellX+0, cellY+1)];
    float velLeft = velocitiesX[idxX(cellX+0, cellY+0)];
    float velRight = velocitiesX[idxX(cellX+1, cellY+0)];
    float velBottom = velocitiesY[idxY(cellX+0, cellY+0)];

    float gradX = (velRight - velLeft) / cellSize;
    float gradY = (velTop - velBottom) / cellSize;

    return gradX + gradY;
}

std::tuple<size_t, size_t> Fluid::findCell(Particle p)
{
    //delicate - explicitly define as tuple
    return {static_cast<size_t>(p.pos.x/cellSize), static_cast<size_t>(p.pos.y/cellSize)};
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

void Fluid::solveIncompressibility(float dt) // Gauss-Seidel
{
    std::fill(pressures.begin(), pressures.end(), 0.0f);

    float cp = density * cellSize / dt;

    for (int iter = 0; iter < numPressureIters; ++iter)
    {
        for (int i = 1; i < size - 1; ++i)
        {
            for (int j = 1; j < size - 1; ++j)
            {
                if (cellType[idxC(i, j)] != FLUID_CELL)
                    continue;

                float sLeft   = s[idxC(i - 1, j)];
                float sRight  = s[idxC(i + 1, j)];
                float sBottom = s[idxC(i, j - 1)];
                float sTop    = s[idxC(i, j + 1)];
                float sSum = sLeft + sRight + sBottom + sTop;

                if (sSum == 0.0f)
                    continue;   //fully boxed in by solids, nothing to solve

                float div = velDivAtCell(i, j);

                float p = -div / sSum;
                p *= overRelaxation;

                pressures[idxC(i, j)] += cp * p;

                //apply correction to the 4 surrounding faces
                velocitiesX[idxX(i, j)]     -= sLeft   * p;
                velocitiesX[idxX(i + 1, j)] += sRight  * p;
                velocitiesY[idxY(i, j)]     -= sBottom * p;
                velocitiesY[idxY(i, j + 1)] += sTop    * p;
            }
        }
    }
}
void Fluid::updateCellType()
{
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
        {
            cellType[idxC(i, j)] = (s[idxC(i, j)] == 0.0f) ? SOLID_CELL : AIR_CELL;
        }

    for (const auto& p : particles)
    {
        auto coords = findCell(p);
        size_t col = std::get<0>(coords);
        size_t row = std::get<1>(coords);
        if (col < static_cast<size_t>(size) && row < static_cast<size_t>(size))
        {
            int idx = idxC(static_cast<int>(col), static_cast<int>(row));
            if (cellType[idx] == AIR_CELL)
            { cellType[idx] = FLUID_CELL; }
        }
    }
}

void Fluid::grid2Particles()
{
    for (auto& p : particles)
    {
        //sample u component
        float gxU = p.pos.x / cellSize;
        float gyU = p.pos.y / cellSize - 0.5f;

        int i0u = clampInt(static_cast<int>(std::floor(gxU)), 0, size);
        int j0u = clampInt(static_cast<int>(std::floor(gyU)), 0, size - 1);
        int i1u = clampInt(i0u + 1, 0, size);
        int j1u = clampInt(j0u + 1, 0, size - 1);

        float txU = std::min(std::max(gxU - std::floor(gxU), 0.0f), 1.0f);
        float tyU = std::min(std::max(gyU - std::floor(gyU), 0.0f), 1.0f);

        float w00u = (1 - txU) * (1 - tyU);
        float w10u = txU * (1 - tyU);
        float w01u = (1 - txU) * tyU;
        float w11u = txU * tyU;

        float picU = velocitiesX[idxX(i0u, j0u)] * w00u + velocitiesX[idxX(i1u, j0u)] * w10u
                     + velocitiesX[idxX(i0u, j1u)] * w01u + velocitiesX[idxX(i1u, j1u)] * w11u;

        float oldU = preVelocitiesX[idxX(i0u, j0u)] * w00u + preVelocitiesX[idxX(i1u, j0u)] * w10u
                     + preVelocitiesX[idxX(i0u, j1u)] * w01u + preVelocitiesX[idxX(i1u, j1u)] * w11u;

        float flipU = p.vel.x + (picU - oldU);

        //sample v component
        float gxV = p.pos.x / cellSize - 0.5f;
        float gyV = p.pos.y / cellSize;

        int i0v = clampInt(static_cast<int>(std::floor(gxV)), 0, size - 1);
        int j0v = clampInt(static_cast<int>(std::floor(gyV)), 0, size);
        int i1v = clampInt(i0v + 1, 0, size - 1);
        int j1v = clampInt(j0v + 1, 0, size);

        float txV = std::min(std::max(gxV - std::floor(gxV), 0.0f), 1.0f);
        float tyV = std::min(std::max(gyV - std::floor(gyV), 0.0f), 1.0f);

        float w00v = (1 - txV) * (1 - tyV);
        float w10v = txV * (1 - tyV);
        float w01v = (1 - txV) * tyV;
        float w11v = txV * tyV;

        float picV = velocitiesY[idxY(i0v, j0v)] * w00v + velocitiesY[idxY(i1v, j0v)] * w10v
                     + velocitiesY[idxY(i0v, j1v)] * w01v + velocitiesY[idxY(i1v, j1v)] * w11v;

        float oldV = preVelocitiesY[idxY(i0v, j0v)] * w00v + preVelocitiesY[idxY(i1v, j0v)] * w10v
                     + preVelocitiesY[idxY(i0v, j1v)] * w01v + preVelocitiesY[idxY(i1v, j1v)] * w11v;

        float flipV = p.vel.y + (picV - oldV);

        //blend
        p.vel.x = flipRatio * flipU + (1.0f - flipRatio) * picU;
        p.vel.y = flipRatio * flipV + (1.0f - flipRatio) * picV;
    }
}

//poorly named step loop
std::vector<int> Fluid::findFluid(float dt)
{
    for (auto& p : particles)
    {
        vec2 gravity{0.0f, 1.5f};
        p.integrate(dt, gravity);
    }

    particles2Grid();

    preVelocitiesX = velocitiesX;
    preVelocitiesY = velocitiesY;

    updateCellType();
    solveIncompressibility(dt);
    grid2Particles();

    return cellType;
}
