// authors: Jan Kaspar (jan.kaspar@gmail.com)

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"
#include "CondFormats/DataRecord/interface/CTPPSOpticsRcd.h"
#include "CondFormats/DataRecord/interface/CTPPSInterpolatedOpticsRcd.h"

#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "CondFormats/CTPPSReadoutObjects/interface/LHCOpticalFunctionsSetCollection.h"
#include "CondFormats/CTPPSReadoutObjects/interface/LHCInterpolatedOpticalFunctionsSetCollection.h"


class CTPPSInterpolatedOpticalFunctionsESSource : public edm::ESProducer
{
  public:
    CTPPSInterpolatedOpticalFunctionsESSource(const edm::ParameterSet&);
    ~CTPPSInterpolatedOpticalFunctionsESSource() override {}

    std::shared_ptr<LHCInterpolatedOpticalFunctionsSetCollection> produce(const CTPPSInterpolatedOpticsRcd&);

  private:
    float currentCrossingAngle_;
    bool currentDataValid_;
    std::shared_ptr<LHCInterpolatedOpticalFunctionsSetCollection> currentData_;
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

CTPPSInterpolatedOpticalFunctionsESSource::CTPPSInterpolatedOpticalFunctionsESSource(const edm::ParameterSet& iConfig) :
  currentCrossingAngle_(-1.),
  currentDataValid_(false)
{
  setWhatProduced(this, &CTPPSInterpolatedOpticalFunctionsESSource::produce);
}

//----------------------------------------------------------------------------------------------------

std::shared_ptr<LHCInterpolatedOpticalFunctionsSetCollection> CTPPSInterpolatedOpticalFunctionsESSource::produce(const CTPPSInterpolatedOpticsRcd& iRecord)
{
  // get the input data
  edm::ESHandle<LHCOpticalFunctionsSetCollection> hOFColl;
  iRecord.getRecord<CTPPSOpticsRcd>().get(hOFColl);

  edm::ESHandle<LHCInfo> hLHCInfo;
  iRecord.getRecord<LHCInfoRcd>().get(hLHCInfo);

  // process
  if (!currentDataValid_ || hLHCInfo->crossingAngle() != currentCrossingAngle_)
  {
    currentCrossingAngle_ = hLHCInfo->crossingAngle();

    if (currentCrossingAngle_ == 0.)
    {
      edm::LogWarning("CTPPSInterpolatedOpticalFunctionsESSource") << "Invalid crossing angle, no optical functions produced.";

      currentDataValid_ = false;
      currentData_ = std::make_shared<LHCInterpolatedOpticalFunctionsSetCollection>();
    } else {
      edm::LogInfo("CTPPSInterpolatedOpticalFunctionsESSource") << "Crossing angle has changed to " << currentCrossingAngle_ << ".";

      if (hOFColl->size() == 1)
      {
        // case with single-xangle input
        const auto &it = hOFColl->begin();
        if (fabs(currentCrossingAngle_ - it->first) < 1e-6)
        {
          currentData_ = std::make_shared<LHCInterpolatedOpticalFunctionsSetCollection>();
          for (const auto &rp_p : it->second)
          {
            const auto rpId = rp_p.first;
            LHCInterpolatedOpticalFunctionsSet iof(rp_p.second);
            iof.initializeSplines();
            currentData_->emplace(rpId, std::move(iof));
          }
        } else {
          throw cms::Exception("CTPPSInterpolatedOpticalFunctionsESSource") << "Cannot interpolate: input given only for xangle "
            << it->first << " while interpolation requested for " << currentCrossingAngle_ << ".";
        }
      } else {
        // case with multi-xangle input

        // find the closest xangle points for interpolation
        auto it1 = hOFColl->begin();
        auto it2 = std::next(it1);

        if (currentCrossingAngle_ > it1->first)
        {
          for (; it1 != hOFColl->end(); ++it1)
          {
            it2 = std::next(it1);

            if (it2 == hOFColl->end())
            {
              it2 = it1;
              it1 = std::prev(it1);
              break;
            }

            if (it1->first <= currentCrossingAngle_ && currentCrossingAngle_ < it2->first)
              break;
          }
        }

        const auto &xangle1 = it1->first;
        const auto &xangle2 = it2->first;

        const auto &ofs1 = it1->second;
        const auto &ofs2 = it2->second;

        // do the interpoaltion RP by RP
        currentData_ = std::make_shared<LHCInterpolatedOpticalFunctionsSetCollection>();
        for (const auto &rp_p : ofs1)
        {
          const auto rpId = rp_p.first;
          const auto &rp_it2 = ofs2.find(rpId);
          if (rp_it2 == ofs2.end())
            throw cms::Exception("CTPPSInterpolatedOpticalFunctionsESSource") << "Mismatch between ofs1 and ofs2.";

          const auto &of1 = rp_p.second;
          const auto &of2 = rp_it2->second;

          LHCInterpolatedOpticalFunctionsSet iof;
          iof.m_z = of1.getScoringPlaneZ();

          const size_t num_xi_vals = of1.getXiValues().size();

          iof.m_xi_values.resize(num_xi_vals);

          for (size_t fi = 0; fi < of1.getFcnValues().size(); ++fi)
          {
            iof.m_fcn_values[fi].resize(num_xi_vals);

            for (size_t pi = 0; pi < num_xi_vals; ++pi)
            {
              double xi = of1.getXiValues()[pi];
              double xi_control = of2.getXiValues()[pi];

              if (fabs(xi - xi_control) > 1e-6)
                throw cms::Exception("CTPPSInterpolatedOpticalFunctionsESSource") << "Xi mismatch between ofs1 and ofs2.";

              iof.m_xi_values[pi] = xi;

              double v1 = of1.getFcnValues()[fi][pi];
              double v2 = of2.getFcnValues()[fi][pi];
              iof.m_fcn_values[fi][pi] = v1 + (v2 - v1) / (xangle2 - xangle1) * (currentCrossingAngle_ - xangle1);
            }
          }

          iof.initializeSplines();

          currentData_->emplace(rpId, std::move(iof));
        }
      }

      currentDataValid_ = true;
    }
  }

  return currentData_;
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_EVENTSETUP_MODULE(CTPPSInterpolatedOpticalFunctionsESSource);
