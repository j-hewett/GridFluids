#pragma once

#include <vector>
#include <tuple>

struct vec2
{
    float x, y;

    vec2(float x = 0.0, float y = 0.0) : x(x), y(y) {}

    vec2 operator+(const vec2& other) const
    {
        return vec2(x + other.x, y + other.y);
    }
    vec2& operator+=(const vec2& other)
    {
        x += other.x;
        y += other.y;
        return *this;
    }

    vec2 operator-(const vec2& other) const
    {
        return vec2(x - other.x, y - other.y);
    }
    vec2 operator-=(const vec2& other)
    {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    vec2 operator*(const float s) const
    {
        return vec2(x * s, y * s);
    }
    vec2 operator*(const int s) const
    {
        return vec2(x * s, y * s);
    }
};

struct Particle
{
    vec2 pos;
    vec2 vel;
    float radius;

    // methods
    Particle() : radius(0.005) {}

    void integrate(float dt, vec2 gravity)
    {
        vel += gravity * dt;
        pos += vel * dt;

        const float minX = 0.0f, maxX = 1.0f;
        const float minY = 0.0f, maxY = 1.0f;

        const float CoR = 0.9f;

        if (pos.x < minX)
        {
            pos.x = minX + radius; vel.x = -vel.x * CoR;
        }
        else if (pos.x > maxX)
        {
            pos.x = maxX - radius; vel.x = -vel.x * CoR;
        }

        if (pos.y < minY)
        {
            pos.y = minY + radius; vel.y = -vel.y * CoR;
        }
        else if (pos.y > maxY) {
            pos.y = maxY - radius; vel.y = -vel.y * CoR;
        }
    }
};

class Fluid
{
public:
    Fluid(int size, int n_particles);
    std::vector<int> simulate(float dt);

    enum CellType { SOLID_CELL = 0, AIR_CELL = 1, FLUID_CELL = 2 };

private:
    void genParticles();
    void initBoundaries();
    std::tuple<size_t, size_t> findCell(Particle p);
    void clearGrid();
    void particles2Grid();
    void updateCellType();
    void solveIncompressibility(float dt);
    void grid2Particles();

    float velDivAtCell(int cellX, int cellY);
    float sampleGrid(const std::vector<float>& grid, float gx, float gy,
                     int iMax, int jMax) const;

    int idxX(int i, int j) const { return i + (size+1)*j; }
    int idxY(int i, int j) const { return i + size*j; }
    int idxC(int i, int j) const { return i + size*j; }

private:
    int size;
    float cellSize;
    std::vector<float> s; //float, used in arithmetic during pressure solve
    std::vector<int> cellType;
    std::vector<Particle> particles;

    std::vector<float> preVelocitiesX, preVelocitiesY;
    std::vector<float> velocitiesX, velocitiesY;
    std::vector<float> weightsX, weightsY;
    std::vector<float> pressures;
    std::vector<float> divergence;

    static constexpr int numPressureIters = 50;
    static constexpr float overRelaxation = 1.9f;
    static constexpr float density = 1000.0f;
    float flipRatio = 0.9f;
};
