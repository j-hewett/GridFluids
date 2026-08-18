#pragma once

#include <vector>
#include <tuple>

class Fluid
{
public:
    Fluid(int size, int n_particles);
    std::vector<int> simulate(float dt);

    enum CellType { SOLID_CELL = 0, AIR_CELL = 1, FLUID_CELL = 2 };

private:
    void genParticles();
    void initBoundaries();
    std::tuple<size_t, size_t> findCell(size_t idx);
    void updateParticleDensity();
    void clearGrid();
    void particles2Grid();
    void pushParticlesApart(int numIters);
    void updateCellType();
    void solveIncompressibility(float dt);
    void grid2Particles();
    void integrateParticles(float dt, float gravity);
    void handleParticleCollisions();

    bool isUFaceValid(int i, int j) const;
    bool isVFaceValid(int i, int j) const;

    float velDivAtCell(int cellX, int cellY);
    float sampleGrid(const std::vector<float>& grid, float gx, float gy,
                     int iMax, int jMax) const;

    int idxX(int i, int j) const { return i + (size+1)*j; }
    int idxY(int i, int j) const { return i + size*j; }
    int idxC(int i, int j) const { return i + size*j; }

private:
    int n_particles;
    int size;
    float cellSize;
    std::vector<float> s; //float, used in arithmetic during pressure solve
    std::vector<int> cellType;
    std::vector<Particle> particles;

    float particleRadius = 0.005f;
    std::vector<float> pPosX;
    std::vector<float> pPosY;
    std::vector<float> pVelX;
    std::vector<float> pVelY;

    std::vector<float> preVelocitiesX, preVelocitiesY;
    std::vector<float> prevScatterVelocitiesX, prevScatterVelocitiesY;
    std::vector<float> velocitiesX, velocitiesY;
    std::vector<float> weightsX, weightsY;
    std::vector<float> pressures;
    std::vector<float> divergence;
    std::vector<float> particleDensity;

    float particleRestDensity = 0.0f;

    //spatial hashing
    float particleSpacing;
    int pNumX, pNumY, pNumCells;
    std::vector<int> numCellParticles;
    std::vector<int> firstCellParticle;
    std::vector<int> cellParticleIds;
    static constexpr int numParticleIters = 2;

    static constexpr int numPressureIters = 50;
    static constexpr float overRelaxation = 1.9f;
    static constexpr float density = 1000.0f;
    static constexpr float driftCompensationK = 1.0f;
    float flipRatio = 0.95f;
};
