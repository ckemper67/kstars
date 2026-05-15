/*
    SPDX-FileCopyrightText: 2024 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include "donutsguider.h"
#include "ekos_guide_debug.h"

namespace Ekos
{

std::vector<double> DonutsGuider::toDoubleBuffer(
    const QSharedPointer<FITSData> &data, double &clip)
{
    const auto   &stats = data->getStatistics();
    const int     n     = data->width() * data->height();
    const uint8_t *raw  = data->getImageBuffer();

    std::vector<double> out(n);

    switch (stats.dataType)
    {
        case TBYTE:
            clip = 255 * 0.95;
            for (int i = 0; i < n; ++i) out[i] = raw[i];
            break;
        case TSHORT: {
            clip = 32767 * 0.95;
            const int16_t *buf = reinterpret_cast<const int16_t *>(raw);
            for (int i = 0; i < n; ++i) out[i] = buf[i];
            break;
        }
        case TUSHORT: {
            clip = 65535 * 0.95;
            const uint16_t *buf = reinterpret_cast<const uint16_t *>(raw);
            for (int i = 0; i < n; ++i) out[i] = buf[i];
            break;
        }
        case TLONG: {
            clip = 2147483647.0 * 0.95;
            const int32_t *buf = reinterpret_cast<const int32_t *>(raw);
            for (int i = 0; i < n; ++i) out[i] = buf[i];
            break;
        }
        case TULONG: {
            clip = 4294967295.0 * 0.95;
            const uint32_t *buf = reinterpret_cast<const uint32_t *>(raw);
            for (int i = 0; i < n; ++i) out[i] = buf[i];
            break;
        }
        case TFLOAT: {
            clip = stats.max[0] * 0.95;
            const float *buf = reinterpret_cast<const float *>(raw);
            for (int i = 0; i < n; ++i) out[i] = buf[i];
            break;
        }
        case TLONGLONG: {
            clip = 9.22e18 * 0.95;
            const int64_t *buf = reinterpret_cast<const int64_t *>(raw);
            for (int i = 0; i < n; ++i) out[i] = buf[i];
            break;
        }
        case TDOUBLE: {
            clip = stats.max[0] * 0.95;
            const double *buf = reinterpret_cast<const double *>(raw);
            for (int i = 0; i < n; ++i) out[i] = buf[i];
            break;
        }
        default:
            clip = stats.max[0] * 0.95;
            qCWarning(KSTARS_EKOS_GUIDE) << "DONUTS: unsupported dataType" << stats.dataType;
            break;
    }
    return out;
}

void DonutsGuider::reset()
{
    QMutexLocker lock(&m_Mutex);
    m_Guider.reset();
}

void DonutsGuider::setReference(const QSharedPointer<FITSData> &data)
{
    QMutexLocker lock(&m_Mutex);
    if (!data) return;
    if (data->getStatistics().width == 0)
        data->calculateStats(true);

    double clip = 0;
    auto pixels = toDoubleBuffer(data, clip);
    m_Guider.setReference(pixels.data(), data->width(), data->height());
    qCDebug(KSTARS_EKOS_GUIDE) << "DONUTS: Reference frame set.";
}

DonutsGuider::Transform DonutsGuider::calculateTransform(
    const QSharedPointer<FITSData> &data)
{
    QMutexLocker lock(&m_Mutex);
    if (!data) return {};
    if (data->getStatistics().width == 0)
        data->calculateStats(true);

    double clip = 0;
    auto pixels = toDoubleBuffer(data, clip);
    return m_Guider.measure(pixels.data(), data->width(), data->height());
}

} // namespace Ekos
