#include "Fluid.h"
#include <random>

static inline int clampInt(int v, int lo, int hi)
{
    return (v < lo) ? lo : (v > hi) ? hi : v;
}

Fluid::Fluid(int size, int n_particles)
    :   n_particles(n_particles),
    size(size),
    cellSize(1.0 / size),
    s(size * size, 1.0f),
    cellType(size * size, AIR_CELL),

    pPosX(n_particles),
    pPosY(n_particles),
    pVelX(n_particles),
    pVelY(n_particles),

    velocitiesX((size+1) * size, 0.0),
    velocitiesY(size * (size+1), 0.0),
    preVelocitiesX((size+1) * size, 0.0f),
    preVelocitiesY(size * (size+1), 0.0f),
    prevScatterVelocitiesX((size+1) * size, 0.0f),
    prevScatterVelocitiesY(size * (size+1), 0.0f),
    pressures(size*size, 0.0),
    particleDensity(size*size, 0.0f),
    weightsX((size+1)*size, 0.0f),
    weightsY(size*(size+1), 0.0f)
{
    genParticles();
    initBoundaries();
    particleRestDensity = 3.0f;

    //init spatial hash relative to particle diameter
    float r = particleRadius;
    particleSpacing = 2.2f * r;
    float pInvSpacing = 1.0f / particleSpacing;
    pNumX = static_cast<int>(1.0f * pInvSpacing) + 1;
    pNumY = static_cast<int>(1.0f * pInvSpacing) + 1;
    pNumCells = pNumX * pNumY;

    numCellParticles.assign(pNumCells, 0);
    firstCellParticle.assign(pNumCells + 1, 0);
    cellParticleIds.assign(n_particles, 0);
}

void Fluid::genParticles()
{
    static std::default_random_engine rng(42);
    std::uniform_real_distribution<float> jitterDist(-0.001f, 0.001f);

    const float targetDensity = std::max(particleRestDensity, 1.0f);
    const float spacing = cellSize / std::sqrt(targetDensity);

    // lay particles out in a roughly square block
    int numX = static_cast<int>(std::sqrt(static_cast<float>(n_particles)));
    if (numX < 1) numX = 1;
    int numY = (static_cast<int>(n_particles) + numX - 1) / numX;

    float blockWidth  = numX * spacing;

    // centered horizontally, flush against the wall nearest pos.y = cellSize
    float startX = std::max(0.5f - blockWidth * 0.5f, cellSize + 0.005f);
    float startY = cellSize + 0.01f;

    size_t idx = 0;
    for (int i = 0; i < numX && idx < n_particles; ++i)
    {
        for (int j = 0; j < numY && idx < n_particles; ++j)
        {
            float jitterX = jitterDist(rng);
            float jitterY = jitterDist(rng);
            pPosX[idx] = startX + i * spacing + jitterX;
            pPosY[idx] = startY + j * spacing + jitterY;
            pVelX[idx] = 0.0f;
            pVelY[idx] = 0.0f;
            ++idx;
        }
    }
}

void Fluid::integrateParticles(float dt, float gravity)
{
    float* __restrict velY = pVelY.data();
    float* __restrict posX = pPosX.data();
    float* __restrict posY = pPosY.data();
    float* __restrict velX = pVelX.data();

    for (int i = 0;i<n_particles; i++)
    {
        velY[i] += gravity * dt;
        posX[i] += velX[i] * dt;
        posY[i] += velY[i] * dt;
    }
}

void Fluid::handleParticleCollisions()
{
    float* __restrict posX = pPosX.data();
    float* __restrict posY = pPosY.data();
    float* __restrict velX = pVelX.data();
    float* __restrict velY = pVelY.data();

    //separate pass so that later hot path can vectorise
    for (int i = 0; i < n_particles; i++)
    {
        bool finite = std::isfinite(posX[i]) && std::isfinite(posY[i]) &&
                      std::isfinite(velX[i]) && std::isfinite(velY[i]);
        if (!finite)
        {
            posX[i] = 0.5f; posY[i] = 0.5f;
            velX[i] = 0.0f; velY[i] = 0.0f;
        }
    }

    const float wallMinX = cellSize, wallMaxX = 1.0f - cellSize;
    const float wallMinY = cellSize, wallMaxY = 1.0f - cellSize;
    const float minX = wallMinX + particleRadius;
    const float maxX = wallMaxX - particleRadius;
    const float minY = wallMinY + particleRadius;
    const float maxY = wallMaxY - particleRadius;

    for (int i = 0; i < n_particles; i++)
    {
        float px = posX[i];
        bool hitX = (px < minX) || (px > maxX);
        float clampedX = (px < minX) ? minX : (px > maxX ? maxX : px);
        posX[i] = clampedX;
        velX[i] = velX[i] * (1.0f - static_cast<float>(hitX));

        float py = posY[i];
        bool hitY = (py < minY) || (py > maxY);
        float clampedY = (py < minY) ? minY : (py > maxY ? maxY : py);
        posY[i] = clampedY;
        velY[i] = velY[i] * (1.0f - static_cast<float>(hitY));
    }
}

void Fluid::pushParticlesApart(int numIters)
{
    float pInvSpacing = 1.0f / particleSpacing;

    std::fill(numCellParticles.begin(), numCellParticles.end(), 0);

    for (int i = 0; i < n_particles; i++)
    {
        int xi = clampInt(static_cast<int>(pPosX[i] * pInvSpacing), 0, pNumX - 1);
        int yi = clampInt(static_cast<int>(pPosY[i] * pInvSpacing), 0, pNumY - 1);
        numCellParticles[xi * pNumY + yi]++;
    }

    int first = 0;
    for (int i = 0; i < pNumCells; ++i)
    {
        first += numCellParticles[i];
        firstCellParticle[i] = first;
    }
    firstCellParticle[pNumCells] = first;

    for (int i = 0; i < n_particles; i++)
    {
        int xi = clampInt(static_cast<int>(pPosX[i] * pInvSpacing), 0, pNumX - 1);
        int yi = clampInt(static_cast<int>(pPosY[i] * pInvSpacing), 0, pNumY - 1);
        int cellNr = xi * pNumY + yi;
        firstCellParticle[cellNr]--;
        cellParticleIds[firstCellParticle[cellNr]] = i;
    }

    float minDist = 2.0f * particleRadius;
    float minDist2 = minDist * minDist;

    for (int iter = 0; iter < numIters; ++iter)
    {
        for (int i = 0; i < n_particles; i++)
        {
            float px = pPosX[i];
            float py = pPosY[i];

            int pxi = static_cast<int>(px * pInvSpacing);
            int pyi = static_cast<int>(py * pInvSpacing);
            int x0 = std::max(pxi - 1, 0);
            int y0 = std::max(pyi - 1, 0);
            int x1 = std::min(pxi + 1, pNumX - 1);
            int y1 = std::min(pyi + 1, pNumY - 1);

            for (int xi = x0; xi <= x1; ++xi)
            {
                for (int yi = y0; yi <= y1; ++yi)
                {
                    int cellNr = xi * pNumY + yi;
                    int first = firstCellParticle[cellNr];
                    int last = firstCellParticle[cellNr + 1];

                    for (int j = first; j < last; ++j)
                    {
                        int id = cellParticleIds[j];
                        if (id == i)
                            continue;

                        float qx = pPosX[id];
                        float qy = pPosY[id];

                        float dx = qx - px;
                        float dy = qy - py;
                        float d2 = dx * dx + dy * dy;
                        if (d2 > minDist2 || d2 == 0.0f)
                            continue;

                        float d = std::sqrt(d2);
                        float s = 0.5f * (minDist - d) / d;
                        dx *= s;
                        dy *= s;

                        pPosX[i]  -= dx;
                        pPosY[i]  -= dy;
                        pPosX[id] += dx;
                        pPosY[id] += dy;
                    }
                }
            }
        }
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

std::tuple<size_t, size_t> Fluid::findCell(size_t idx)
{
    int col = clampInt(static_cast<int>(pPosX[idx] / cellSize), 0, size - 1);
    int row = clampInt(static_cast<int>(pPosY[idx] / cellSize), 0, size - 1);
    return {static_cast<size_t>(col), static_cast<size_t>(row)};
}

void Fluid::updateParticleDensity()
{
    std::fill(particleDensity.begin(), particleDensity.end(), 0.0f);

    const float* __restrict posX = pPosX.data();
    const float* __restrict posY = pPosY.data();
    float* __restrict pDens = particleDensity.data();

    float h = cellSize;
    float h2 = 0.5f * h;

    for (int i = 0; i < n_particles; i++)
    {
        float x = std::min(std::max(posX[i], h), (size - 1) * h);
        float y = std::min(std::max(posY[i], h), (size - 1) * h);

        int x0 = static_cast<int>((x - h2) / h);
        float tx = ((x - h2) - x0 * h) / h;
        int x1 = std::min(x0 + 1, size - 2);

        int y0 = static_cast<int>((y - h2) / h);
        float ty = ((y - h2) - y0 * h) / h;
        int y1 = std::min(y0 + 1, size - 2);

        float sx = 1.0f - tx;
        float sy = 1.0f - ty;

        if (x0 < size && y0 < size) pDens[idxC(x0, y0)] += sx * sy;
        if (x1 < size && y0 < size) pDens[idxC(x1, y0)] += tx * sy;
        if (x1 < size && y1 < size) pDens[idxC(x1, y1)] += tx * ty;
        if (x0 < size && y1 < size) pDens[idxC(x0, y1)] += sx * ty;
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

    // u component
    for (int i = 0; i < n_particles; i++)
    {
        float gx = pPosX[i] / cellSize;
        float gy = pPosY[i] / cellSize - 0.5f;

        int i0 = static_cast<int>(std::floor(gx));
        int j0 = static_cast<int>(std::floor(gy));
        int i1 = i0 + 1;
        int j1 = j0 + 1;

        float tx = gx - i0;
        float ty = gy - j0;

        i0 = clampInt(i0, 0, size);
        i1 = clampInt(i1, 0, size);
        j0 = clampInt(j0, 0, size - 1);
        j1 = clampInt(j1, 0, size - 1);

        float w00 = (1 - tx) * (1 - ty);
        float w10 = tx * (1 - ty);
        float w01 = (1 - tx) * ty;
        float w11 = tx * ty;

        float vx = pVelX[i];
        velocitiesX[idxX(i0, j0)] += vx * w00;  weightsX[idxX(i0, j0)] += w00;
        velocitiesX[idxX(i1, j0)] += vx * w10;  weightsX[idxX(i1, j0)] += w10;
        velocitiesX[idxX(i0, j1)] += vx * w01;  weightsX[idxX(i0, j1)] += w01;
        velocitiesX[idxX(i1, j1)] += vx * w11;  weightsX[idxX(i1, j1)] += w11;
    }

    // v component
    for (int i = 0; i < n_particles; i++)
    {
        float gx = pPosX[i] / cellSize - 0.5f;
        float gy = pPosY[i] / cellSize;

        int i0 = static_cast<int>(std::floor(gx));
        int j0 = static_cast<int>(std::floor(gy));
        int i1 = i0 + 1;
        int j1 = j0 + 1;

        float tx = gx - i0;
        float ty = gy - j0;

        i0 = clampInt(i0, 0, size - 1);
        i1 = clampInt(i1, 0, size - 1);
        j0 = clampInt(j0, 0, size);
        j1 = clampInt(j1, 0, size);

        float w00 = (1 - tx) * (1 - ty);
        float w10 = tx * (1 - ty);
        float w01 = (1 - tx) * ty;
        float w11 = tx * ty;

        float vy = pVelY[i];
        velocitiesY[idxY(i0, j0)] += vy * w00;  weightsY[idxY(i0, j0)] += w00;
        velocitiesY[idxY(i1, j0)] += vy * w10;  weightsY[idxY(i1, j0)] += w10;
        velocitiesY[idxY(i0, j1)] += vy * w01;  weightsY[idxY(i0, j1)] += w01;
        velocitiesY[idxY(i1, j1)] += vy * w11;  weightsY[idxY(i1, j1)] += w11;
    }

    for (size_t idx = 0; idx < velocitiesX.size(); ++idx)
        if (weightsX[idx] > 0.0f)
            velocitiesX[idx] /= weightsX[idx];

    for (size_t idx = 0; idx < velocitiesY.size(); ++idx)
        if (weightsY[idx] > 0.0f)
            velocitiesY[idx] /= weightsY[idx];

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

void Fluid::solveIncompressibility(float dt)
{
    std::fill(pressures.begin(), pressures.end(), 0.0f);

    preVelocitiesX = velocitiesX;
    preVelocitiesY = velocitiesY;

    float cp = density * cellSize / dt;

    // color: 0 = red (i+j even), 1 = black (i+j odd)
    for (int iter = 0; iter < numPressureIters; ++iter)
    {
        for (int color = 0; color < 2; ++color)
        {
            for (int i = 1; i < size - 1; ++i)
            {
                // start j so that (i+j) matches this color, then step by 2
                int jStart = 1 + ((i + 1 + color) % 2);

                for (int j = jStart; j < size - 1; j += 2)
                {
                    if (cellType[idxC(i, j)] != FLUID_CELL)
                        continue;

                    float sLeft   = s[idxC(i - 1, j)];
                    float sRight  = s[idxC(i + 1, j)];
                    float sBottom = s[idxC(i, j - 1)];
                    float sTop    = s[idxC(i, j + 1)];
                    float sSum = sLeft + sRight + sBottom + sTop;

                    if (sSum == 0.0f)
                        continue;

                    float div = velDivAtCell(i, j);

                    if (particleRestDensity > 0.0f)
                    {
                        float compression = particleDensity[idxC(i, j)] - particleRestDensity;
                        if (compression > 0.0f)
                            div -= driftCompensationK * compression;
                    }

                    float p = -div / sSum;
                    p *= overRelaxation;

                    pressures[idxC(i, j)] += cp * p;

                    velocitiesX[idxX(i, j)]     -= sLeft   * p;
                    velocitiesX[idxX(i + 1, j)] += sRight  * p;
                    velocitiesY[idxY(i, j)]     -= sBottom * p;
                    velocitiesY[idxY(i, j + 1)] += sTop    * p;
                }
            }
        }
    }

    const float maxVel = 10.0f * cellSize / dt;
    for (auto& v : velocitiesX) v = std::clamp(v, -maxVel, maxVel);
    for (auto& v : velocitiesY) v = std::clamp(v, -maxVel, maxVel);
}

void Fluid::updateCellType()
{
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            cellType[idxC(i, j)] = (s[idxC(i, j)] == 0.0f) ? SOLID_CELL : AIR_CELL;

    for (int p = 0; p < n_particles; p++)
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
    for (int i = 0; i < n_particles; i++)
    {
        // sample u component
        float gxU = pPosX[i] / cellSize;
        float gyU = pPosY[i] / cellSize - 0.5f;

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
        float picU = pVelX[i];
        if (wSumU > 0.0f)
        {
            picU = (v00*w00u*velocitiesX[idxX(i0u,j0u)] + v10*w10u*velocitiesX[idxX(i1u,j0u)]
                    + v01*w01u*velocitiesX[idxX(i0u,j1u)] + v11*w11u*velocitiesX[idxX(i1u,j1u)]) / wSumU;

            float oldU = (v00*w00u*preVelocitiesX[idxX(i0u,j0u)] + v10*w10u*preVelocitiesX[idxX(i1u,j0u)]
                          + v01*w01u*preVelocitiesX[idxX(i0u,j1u)] + v11*w11u*preVelocitiesX[idxX(i1u,j1u)]) / wSumU;

            flipU = pVelX[i] + (picU - oldU);
        }
        else
        {
            flipU = pVelX[i];
        }

        // sample v component
        float gxV = pPosX[i] / cellSize - 0.5f;
        float gyV = pPosY[i] / cellSize;

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
        float picV = pVelY[i];
        if (wSumV > 0.0f)
        {
            picV = (u00*w00v*velocitiesY[idxY(i0v,j0v)] + u10*w10v*velocitiesY[idxY(i1v,j0v)]
                    + u01*w01v*velocitiesY[idxY(i0v,j1v)] + u11*w11v*velocitiesY[idxY(i1v,j1v)]) / wSumV;

            float oldV = (u00*w00v*preVelocitiesY[idxY(i0v,j0v)] + u10*w10v*preVelocitiesY[idxY(i1v,j0v)]
                          + u01*w01v*preVelocitiesY[idxY(i0v,j1v)] + u11*w11v*preVelocitiesY[idxY(i1v,j1v)]) / wSumV;

            flipV = pVelY[i] + (picV - oldV);
        }
        else
        {
            flipV = pVelY[i];
        }

        pVelX[i] = flipRatio * flipU + (1.0f - flipRatio) * picU;
        pVelY[i] = flipRatio * flipV + (1.0f - flipRatio) * picV;
    }
}

std::vector<int> Fluid::simulate(float dt)
{
    float gravity = 1.0f;
    integrateParticles(dt, gravity);
    pushParticlesApart(numParticleIters);
    handleParticleCollisions();

    particles2Grid();
    updateCellType();
    updateParticleDensity();
    solveIncompressibility(dt);
    grid2Particles();

    return cellType;
}