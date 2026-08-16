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
    prevScatterVelocitiesX((size+1) * size, 0.0f),
    prevScatterVelocitiesY(size * (size+1), 0.0f),
    pressures(size*size, 0.0),
    weightsX((size+1)*size, 0.0f),
    weightsY(size*(size+1), 0.0f)
{
    genParticles();
    initBoundaries();
}


void Fluid::genParticles()
{
    static std::default_random_engine rng(42);
    std::uniform_real_distribution<float> dist(cellSize + 0.005f, 1.0f - cellSize - 0.005f);

    for (auto& p : particles)
    {
        p.pos = vec2(dist(rng), dist(rng));
        vec2 randVel = vec2((dist(rng)-0.5)*0.5, (dist(rng)-0.5)*0.5);
        p.vel = randVel;
    }
}

void Fluid::integrateParticles(float dt, vec2 gravity)
{
    for (auto& p : particles)
    {
        p.vel += gravity * dt;
        p.pos += p.vel * dt;
    }
}

void Fluid::handleParticleCollisions()
{
    const float wallMinX = cellSize, wallMaxX = 1.0f - cellSize;
    const float wallMinY = cellSize, wallMaxY = 1.0f - cellSize;

    for (auto& p : particles)
    {
        if (!std::isfinite(p.pos.x) || !std::isfinite(p.pos.y) ||
            !std::isfinite(p.vel.x) || !std::isfinite(p.vel.y))
        {
            p.pos = vec2(0.5f, 0.5f);
            p.vel = vec2(0.0f, 0.0f);
        }
        const float minX = wallMinX + p.radius;
        const float maxX = wallMaxX - p.radius;
        const float minY = wallMinY + p.radius;
        const float maxY = wallMaxY - p.radius;

        if (p.pos.x < minX)      { p.pos.x = minX; p.vel.x = 0; }
        else if (p.pos.x > maxX) { p.pos.x = maxX; p.vel.x = 0; }

        if (p.pos.y < minY)      { p.pos.y = minY; p.vel.y = 0; }
        else if (p.pos.y > maxY) { p.pos.y = maxY; p.vel.y = 0; }
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

//returns true if this u-face (i,j) is adjacent to at least one non-air cell
bool Fluid::isUFaceValid(int i, int j) const
{
    bool valid = false;
    if (i <= size - 1) valid = valid || (cellType[idxC(i, j)] != AIR_CELL);
    if (i >= 1)         valid = valid || (cellType[idxC(i - 1, j)] != AIR_CELL);
    return valid;
}

//returns true if this v-face (i,j) is adjacent to at least one non-air cell
bool Fluid::isVFaceValid(int i, int j) const
{
    bool valid = false;
    if (j <= size - 1) valid = valid || (cellType[idxC(i, j)] != AIR_CELL);
    if (j >= 1)         valid = valid || (cellType[idxC(i, j - 1)] != AIR_CELL);
    return valid;
}

float Fluid::velDivAtCell(int cellX, int cellY)
{
    if (cellType[idxC(cellX, cellY)] != FLUID_CELL)
        return 0.0f;

    float velTop = velocitiesY[idxY(cellX+0, cellY+1)];
    float velLeft = velocitiesX[idxX(cellX+0, cellY+0)];
    float velRight = velocitiesX[idxX(cellX+1, cellY+0)];
    float velBottom = velocitiesY[idxY(cellX+0, cellY+0)];

    float gradX = (velRight - velLeft);
    float gradY = (velTop - velBottom);

    return gradX + gradY;
}

std::tuple<size_t, size_t> Fluid::findCell(Particle p)
{
    int col = clampInt(static_cast<int>(p.pos.x / cellSize), 0, size - 1);
    int row = clampInt(static_cast<int>(p.pos.y / cellSize), 0, size - 1);
    return {static_cast<size_t>(col), static_cast<size_t>(row)};
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

    //restore velocity on faces touching solid cells
    for (int i = 0; i <= size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            bool solidRight = (i <= size - 1) ? (s[idxC(i, j)] == 0.0f) : true;
            bool solidLeft  = (i >= 1)        ? (s[idxC(i - 1, j)] == 0.0f) : true;
            if (solidLeft || solidRight)
                velocitiesX[idxX(i, j)] = prevScatterVelocitiesX[idxX(i, j)];
        }
    }

    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j <= size; ++j)
        {
            bool solidTop    = (j <= size - 1) ? (s[idxC(i, j)] == 0.0f) : true;
            bool solidBottom = (j >= 1)        ? (s[idxC(i, j - 1)] == 0.0f) : true;
            if (solidTop || solidBottom)
                velocitiesY[idxY(i, j)] = prevScatterVelocitiesY[idxY(i, j)];
        }
    }
}

void Fluid::solveIncompressibility(float dt) // Gauss-Seidel
{
    std::fill(pressures.begin(), pressures.end(), 0.0f);

    preVelocitiesX = velocitiesX;
    preVelocitiesY = velocitiesY;

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

    const float maxVel = 10.0f * cellSize / dt; // ~10 cells/frame ceiling
    for (auto& v : velocitiesX) v = std::clamp(v, -maxVel, maxVel);
    for (auto& v : velocitiesY) v = std::clamp(v, -maxVel, maxVel);
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

        float v00 = isUFaceValid(i0u, j0u) ? 1.0f : 0.0f;
        float v10 = isUFaceValid(i1u, j0u) ? 1.0f : 0.0f;
        float v01 = isUFaceValid(i0u, j1u) ? 1.0f : 0.0f;
        float v11 = isUFaceValid(i1u, j1u) ? 1.0f : 0.0f;

        float wSumU = v00*w00u + v10*w10u + v01*w01u + v11*w11u;

        float flipU;
        float picU = p.vel.x;
        if (wSumU > 0.0f)
        {
            picU = (v00*w00u*velocitiesX[idxX(i0u,j0u)] + v10*w10u*velocitiesX[idxX(i1u,j0u)]
                          + v01*w01u*velocitiesX[idxX(i0u,j1u)] + v11*w11u*velocitiesX[idxX(i1u,j1u)]) / wSumU;

            float oldU = (v00*w00u*preVelocitiesX[idxX(i0u,j0u)] + v10*w10u*preVelocitiesX[idxX(i1u,j0u)]
                          + v01*w01u*preVelocitiesX[idxX(i0u,j1u)] + v11*w11u*preVelocitiesX[idxX(i1u,j1u)]) / wSumU;

            flipU = p.vel.x + (picU - oldU);
        }
        else
        {
            flipU = p.vel.x;
        }

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

        float u00 = isVFaceValid(i0v, j0v) ? 1.0f : 0.0f;
        float u10 = isVFaceValid(i1v, j0v) ? 1.0f : 0.0f;
        float u01 = isVFaceValid(i0v, j1v) ? 1.0f : 0.0f;
        float u11 = isVFaceValid(i1v, j1v) ? 1.0f : 0.0f;

        float wSumV = u00*w00v + u10*w10v + u01*w01v + u11*w11v;

        float flipV;
        float picV = p.vel.y;
        if (wSumV > 0.0f)
        {
            picV = (u00*w00v*velocitiesY[idxY(i0v,j0v)] + u10*w10v*velocitiesY[idxY(i1v,j0v)]
                          + u01*w01v*velocitiesY[idxY(i0v,j1v)] + u11*w11v*velocitiesY[idxY(i1v,j1v)]) / wSumV;

            float oldV = (u00*w00v*preVelocitiesY[idxY(i0v,j0v)] + u10*w10v*preVelocitiesY[idxY(i1v,j0v)]
                          + u01*w01v*preVelocitiesY[idxY(i0v,j1v)] + u11*w11v*preVelocitiesY[idxY(i1v,j1v)]) / wSumV;

            flipV = p.vel.y + (picV - oldV);
        }
        else
        {
            flipV = p.vel.y;
        }

        //blend
        p.vel.x = flipRatio * flipU + (1.0f - flipRatio) * picU;
        p.vel.y = flipRatio * flipV + (1.0f - flipRatio) * picV;
    }
}

std::vector<int> Fluid::simulate(float dt)
{
    vec2 gravity = {0.0, -1.5f};
    integrateParticles(dt, gravity);
    handleParticleCollisions();

    particles2Grid();
    updateCellType();
    solveIncompressibility(dt);
    grid2Particles();

    return cellType;
}
