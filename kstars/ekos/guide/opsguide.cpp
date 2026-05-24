/*
    SPDX-FileCopyrightText: 2003 Jasem Mutlaq <mutlaqja@ikarustech.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include "opsguide.h"

#include "Options.h"
#include "kstars.h"
#include "auxiliary/ksnotification.h"
#include "internalguide/internalguider.h"
#include "ekos/auxiliary/stellarsolverprofileeditor.h"
#include "kspaths.h"

#include <KConfigDialog>

#include <QCheckBox>
#include <QComboBox>
#include <QFileDialog>
#include <QPushButton>
#include <QStringList>

namespace Ekos
{
OpsGuide::OpsGuide() : QFrame(KStars::Instance())
{
    setupUi(this);

    //Get a pointer to the KConfigDialog
    m_ConfigDialog = KConfigDialog::exists("guidesettings");

    editGuideProfile->setIcon(QIcon::fromTheme("document-edit"));
    editGuideProfile->setAttribute(Qt::WA_LayoutUsesWidgetRect);

    if (Options::gPGEnabled())
    {
        // this option is an old one. Allowing users who have this set to keep their setting.
        Options::setGPGEnabled(false);
        Options::setRAGuidePulseAlgorithm(GPG_ALGORITHM);
        kcfg_RAGuidePulseAlgorithm->setCurrentIndex(static_cast<int>(GPG_ALGORITHM));
    }

    // Connect comboboxes for algorithm changes
    connect(kcfg_RAGuidePulseAlgorithm, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &OpsGuide::setRAGuidePulseAlg);
    connect(kcfg_DECGuidePulseAlgorithm, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &OpsGuide::setDECGuidePulseAlg);

    // Connect Harmonic Profile checkboxes
    connect(kcfg_RAMPCHarmonicProfile, &QCheckBox::toggled, this, &OpsGuide::slotRAMPCHarmonicToggled);
    connect(kcfg_DECMPCHarmonicProfile, &QCheckBox::toggled, this, &OpsGuide::slotDECMPCHarmonicToggled);

    connect(editGuideProfile, &QAbstractButton::clicked, this, [this]()
    {
        KConfigDialog *optionsEditor = new KConfigDialog(this, "OptionsProfileEditor", Options::self());
        optionsProfileEditor = new StellarSolverProfileEditor(this, Ekos::GuideProfiles, optionsEditor);
#ifdef Q_OS_MACOS
        optionsEditor->setWindowFlags(Qt::Tool | Qt::WindowStaysOnTopHint);
#endif
        KPageWidgetItem *mainPage = optionsEditor->addPage(optionsProfileEditor, i18n("Guide Options Profile Editor"));
        mainPage->setIcon(QIcon::fromTheme("configure"));
        connect(optionsProfileEditor, &StellarSolverProfileEditor::optionsProfilesUpdated, this, &OpsGuide::loadOptionsProfiles);
        optionsProfileEditor->loadProfile(kcfg_GuideOptionsProfile->currentText());
        optionsEditor->show();
    });

    loadOptionsProfiles();

    // Initialize UI states based on active algorithms and profiles
    setRAGuidePulseAlg(kcfg_RAGuidePulseAlgorithm->currentIndex());
    setDECGuidePulseAlg(kcfg_DECGuidePulseAlgorithm->currentIndex());
    slotRAMPCHarmonicToggled(kcfg_RAMPCHarmonicProfile->isChecked());
    slotDECMPCHarmonicToggled(kcfg_DECMPCHarmonicProfile->isChecked());

    connect(m_ConfigDialog, &KConfigDialog::settingsChanged, this, &OpsGuide::settingsUpdated);
}

void OpsGuide::setRAGuidePulseAlg(int index)
{
    switch (index)
    {
        case STANDARD_ALGORITHM:
            kcfg_RAHysteresis->setDisabled(true);
            kcfg_RAIntegralGain->setDisabled(false);
            kcfg_RAProportionalGain->setEnabled(true);
            kcfg_RAMPCTrackingWeight->setDisabled(true);
            kcfg_RAMPCControlPenalty->setDisabled(true);
            kcfg_RAMPCHarmonicProfile->setDisabled(true);
            break;
        case HYSTERESIS_ALGORITHM:
            kcfg_RAHysteresis->setDisabled(false);
            kcfg_RAIntegralGain->setDisabled(true);
            kcfg_RAProportionalGain->setEnabled(true);
            kcfg_RAMPCTrackingWeight->setDisabled(true);
            kcfg_RAMPCControlPenalty->setDisabled(true);
            kcfg_RAMPCHarmonicProfile->setDisabled(true);
            break;
        case LINEAR_ALGORITHM:
        case GPG_ALGORITHM:
            kcfg_RAHysteresis->setDisabled(true);
            kcfg_RAIntegralGain->setDisabled(true);
            kcfg_RAProportionalGain->setEnabled(true);
            kcfg_RAMPCTrackingWeight->setDisabled(true);
            kcfg_RAMPCControlPenalty->setDisabled(true);
            kcfg_RAMPCHarmonicProfile->setDisabled(true);
            break;
        case MPC_ALGORITHM:
            kcfg_RAHysteresis->setDisabled(true);
            kcfg_RAIntegralGain->setDisabled(true);
            kcfg_RAProportionalGain->setDisabled(true);
            kcfg_RAMPCTrackingWeight->setEnabled(!kcfg_RAMPCHarmonicProfile->isChecked());
            kcfg_RAMPCControlPenalty->setEnabled(!kcfg_RAMPCHarmonicProfile->isChecked());
            kcfg_RAMPCHarmonicProfile->setEnabled(true);
            break;
        default:
            kcfg_RAHysteresis->setDisabled(false);
            kcfg_RAIntegralGain->setDisabled(false);
            kcfg_RAProportionalGain->setEnabled(true);
            kcfg_RAMPCTrackingWeight->setDisabled(true);
            kcfg_RAMPCControlPenalty->setDisabled(true);
            kcfg_RAMPCHarmonicProfile->setDisabled(true);
            break;
    }
}

void OpsGuide::setDECGuidePulseAlg(int index)
{
    switch (index)
    {
        case STANDARD_ALGORITHM:
            kcfg_DECHysteresis->setDisabled(true);
            kcfg_DECIntegralGain->setDisabled(false);
            kcfg_DECProportionalGain->setEnabled(true);
            kcfg_DECMPCTrackingWeight->setDisabled(true);
            kcfg_DECMPCControlPenalty->setDisabled(true);
            kcfg_DECMPCHarmonicProfile->setDisabled(true);
            break;
        case HYSTERESIS_ALGORITHM:
            kcfg_DECHysteresis->setDisabled(false);
            kcfg_DECIntegralGain->setDisabled(true);
            kcfg_DECProportionalGain->setEnabled(true);
            kcfg_DECMPCTrackingWeight->setDisabled(true);
            kcfg_DECMPCControlPenalty->setDisabled(true);
            kcfg_DECMPCHarmonicProfile->setDisabled(true);
            break;
        case LINEAR_ALGORITHM:
            kcfg_DECHysteresis->setDisabled(true);
            kcfg_DECIntegralGain->setDisabled(true);
            kcfg_DECProportionalGain->setEnabled(true);
            kcfg_DECMPCTrackingWeight->setDisabled(true);
            kcfg_DECMPCControlPenalty->setDisabled(true);
            kcfg_DECMPCHarmonicProfile->setDisabled(true);
            break;
        case MPC_ALGORITHM:
            kcfg_DECHysteresis->setDisabled(true);
            kcfg_DECIntegralGain->setDisabled(true);
            kcfg_DECProportionalGain->setDisabled(true);
            kcfg_DECMPCTrackingWeight->setEnabled(!kcfg_DECMPCHarmonicProfile->isChecked());
            kcfg_DECMPCControlPenalty->setEnabled(!kcfg_DECMPCHarmonicProfile->isChecked());
            kcfg_DECMPCHarmonicProfile->setEnabled(true);
            break;
        default:
            kcfg_DECHysteresis->setDisabled(false);
            kcfg_DECIntegralGain->setDisabled(false);
            kcfg_DECProportionalGain->setEnabled(true);
            kcfg_DECMPCTrackingWeight->setDisabled(true);
            kcfg_DECMPCControlPenalty->setDisabled(true);
            kcfg_DECMPCHarmonicProfile->setDisabled(true);
            break;
    }
}

void OpsGuide::slotRAMPCHarmonicToggled(bool checked)
{
    if (checked)
    {
        kcfg_RAMPCTrackingWeight->setValue(10.0);
        kcfg_RAMPCControlPenalty->setValue(1.0);
    }
    // Enable/disable spinboxes based on checkbox state and algorithm
    bool isMPC = (kcfg_RAGuidePulseAlgorithm->currentIndex() == MPC_ALGORITHM);
    kcfg_RAMPCTrackingWeight->setEnabled(isMPC && !checked);
    kcfg_RAMPCControlPenalty->setEnabled(isMPC && !checked);
}

void OpsGuide::slotDECMPCHarmonicToggled(bool checked)
{
    if (checked)
    {
        kcfg_DECMPCTrackingWeight->setValue(10.0);
        kcfg_DECMPCControlPenalty->setValue(1.0);
    }
    // Enable/disable spinboxes based on checkbox state and algorithm
    bool isMPC = (kcfg_DECGuidePulseAlgorithm->currentIndex() == MPC_ALGORITHM);
    kcfg_DECMPCTrackingWeight->setEnabled(isMPC && !checked);
    kcfg_DECMPCControlPenalty->setEnabled(isMPC && !checked);
}

void OpsGuide::loadOptionsProfiles()
{
    QString savedOptionsProfiles = QDir(KSPaths::writableLocation(
                                            QStandardPaths::AppLocalDataLocation)).filePath("SavedGuideProfiles.ini");
    if(QFile(savedOptionsProfiles).exists())
        optionsList = StellarSolver::loadSavedOptionsProfiles(savedOptionsProfiles);
    else
        optionsList = getDefaultGuideOptionsProfiles();
    kcfg_GuideOptionsProfile->clear();
    for(SSolver::Parameters param : optionsList)
        kcfg_GuideOptionsProfile->addItem(param.listName);
    kcfg_GuideOptionsProfile->setCurrentIndex(Options::guideOptionsProfile());
}
}
