#pragma once

#include <QWidget>
#include <QPainter>
#include <QPaintEvent>
#include <QDebug>

#include <vector>

class GridWidget : public QWidget
{
    Q_OBJECT

public:
    explicit GridWidget(int rows, int cols, QWidget *parent = nullptr)
        : QWidget(parent)
        , m_rows(rows)
        , m_cols(cols)
        , m_cells(rows * cols, 0)
    {
        if (rows <= 0 || cols <= 0)
        {
            qCritical() << "[GridWidget] Invalid grid dimensions:"
                        << rows << "x" << cols
                        << "- both must be > 0";
        }
        setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    }

    GridWidget(int windowW, int windowH, int cellSizePx, QWidget *parent = nullptr)
        : QWidget(parent)
        , m_rows(windowH / cellSizePx)
        , m_cols(windowW / cellSizePx)
        , m_cells((windowH / cellSizePx) * (windowW / cellSizePx), 0)
    {
        if (cellSizePx <= 0)
        {
            qCritical() << "[GridWidget] Invalid cell size:" << cellSizePx
                        << "- must be > 0";
        }
        if (windowW <= 0 || windowH <= 0)
        {
            qCritical() << "[GridWidget] Invalid window size:"
                        << windowW << "x" << windowH
                        << "- both must be > 0";
        }
        if (m_rows <= 0 || m_cols <= 0)
        {
            qCritical() << "[GridWidget] Calculated grid dimensions:"
                        << m_rows << "x" << m_cols
                        << "- cell size too large for window";
        }

        qInfo() << "[GridWidget] Window" << windowW << "x" << windowH
                << "cellSize" << cellSizePx << "px -> grid"
                << m_cols << "x" << m_rows;

        resize(windowW, windowH);
        setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    }

    void setCells(const std::vector<int>& cells)
    {
        const size_t expected = static_cast<size_t>(m_rows) * m_cols;
        if (cells.size() != expected)
        {
            qWarning() << "[GridWidget] setCells: expected" << expected
                       << "values (" << m_rows << "x" << m_cols << ") but got"
                       << cells.size() << "- ignoring update";
            return;
        }
        m_cells = cells;
        update();
    }

    void setCell(int row, int col, int value)
    {
        if (row < 0 || row >= m_rows || col < 0 || col >= m_cols)
        {
            qWarning() << "[GridWidget] setCell: (" << row << "," << col
                       << ") out of bounds for" << m_rows << "x" << m_cols << "grid";
            return;
        }
        m_cells[row * m_cols + col] = value;
    }

    void setGridSize(int rows, int cols)
    {
        if (rows <= 0 || cols <= 0)
        {
            qCritical() << "[GridWidget] setGridSize: invalid dimensions"
                        << rows << "x" << cols;
            return;
        }
        m_rows = rows;
        m_cols = cols;
        m_cells.assign(rows * cols, 0);
        update();
    }

    int rows() const { return m_rows; }
    int cols() const { return m_cols; }

protected:
    void paintEvent(QPaintEvent * /*event*/) override
    {
        QPainter painter(this);
        painter.setRenderHint(QPainter::Antialiasing, false);

        const float cellSize = qMin(
            static_cast<float>(width())  / m_cols,
            static_cast<float>(height()) / m_rows
            );

        const float gridW = cellSize * m_cols;
        const float gridH = cellSize * m_rows;

        const float offsetX = (width()  - gridW) * 0.5f;
        const float offsetY = (height() - gridH) * 0.5f;

        painter.fillRect(rect(), palette().window());
        painter.save();
        painter.translate(offsetX, offsetY);

        const QColor fillColor(63, 216, 216);
        painter.setPen(Qt::NoPen);
        painter.setBrush(fillColor);

        for (int r = 0; r < m_rows; ++r)
        {
            for (int c = 0; c < m_cols; ++c)
            {
                if (m_cells[r * m_cols + c])
                {
                    painter.drawRect(QRectF(c * cellSize, r * cellSize, cellSize, cellSize));
                }
            }
        }

        QPen gridPen(QColor(200, 200, 200), 0.1, Qt::SolidLine);
        painter.setPen(gridPen);
        painter.setBrush(Qt::NoBrush);

        for (int c = 0; c <= m_cols; ++c)
        {
            const float x = c * cellSize;
            painter.drawLine(QPointF(x, 0), QPointF(x, gridH));
        }

        for (int r = 0; r <= m_rows; ++r)
        {
            const float y = r * cellSize;
            painter.drawLine(QPointF(0, y), QPointF(gridW, y));
        }

        painter.restore();
        painter.end();
    }

private:
    int m_rows;
    int m_cols;
    std::vector<int> m_cells;
};
