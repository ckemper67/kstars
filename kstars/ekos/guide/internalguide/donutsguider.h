/*
    SPDX-FileCopyrightText: 2024 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once

#include <QSharedPointer>
#include <QMutex>
#include "fitsviewer/fitsdata.h"
#include "ekos/guide/donuts/donuts.h"

namespace Ekos
{

/**
 * @class DonutsGuider
 * @short Qt/FITSData wrapper around the Donuts guiding library.
 *
 * Converts FITSData frames to raw double buffers and delegates all
 * algorithm work to Donuts::Guider (see ekos/guide/donuts/donuts.h).
 */
class DonutsGuider
{
    public:
        using Transform = Donuts::Transform;

        DonutsGuider()  = default;
        ~DonutsGuider() = default;

        void setReference(const QSharedPointer<FITSData> &data);
        Transform calculateTransform(const QSharedPointer<FITSData> &data);
        void reset();
        bool hasReference() const { return m_Guider.hasReference(); }

    private:
        static std::vector<double> toDoubleBuffer(
            const QSharedPointer<FITSData> &data, double &clip);

        Donuts::Guider m_Guider;
        QMutex         m_Mutex;

        Q_DISABLE_COPY(DonutsGuider)
};

} // namespace Ekos
