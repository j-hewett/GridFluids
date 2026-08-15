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
private:
    int size;
    float cellSize;
    std::vector<int> gridState;
    std::vector<Particle> particles;

    std::vector<float> velocitiesX, velocitiesY;
    std::vector<float> weightsX, weightsY;
    std::vector<float> pressures;
    std::vector<float> divergence;


public:
    Fluid(int size, int n_particles);

    float velDivAtCell(int cellX, int cellY);
    std::vector<int> findFluid(float dt);

private:
    void genParticles();
    std::tuple<size_t, size_t> findCell(Particle p);
    void clearGrid();
    void particles2Grid();


private:
    int idxX(int i, int j) const { return i + (size+1)*j; }
    int idxY(int i, int j) const { return i + size*j; }
    int idxC(int i, int j) const { return i + size*j; }

};
