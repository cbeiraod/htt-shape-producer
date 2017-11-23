# -*- coding: utf-8 -*-

import copy
import os

from estimation_methods import EstimationMethod, SStoOSEstimationMethod, ABCDEstimationMethod
from histogram import *
from cutstring import *
from systematics import *
from systematic_variations import *
from era import log_query

import logging
logger = logging.getLogger(__name__)


class DataEstimation(EstimationMethod):
    def __init__(self, era, directory, channel):
        super(DataEstimation, self).__init__(
            name="data_obs",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign=None)
        self._channel = channel

    def get_files(self):
        return self.artus_file_names(self.era.data_files(self._channel))

    def get_cuts(self):
        return Cuts()


class HTTEstimation(EstimationMethod):
    def __init__(self, era, directory, channel):
        super(HTTEstimation, self).__init__(
            name="HTT",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_weights(self):
        return Weights(
            Weight("eventWeight", "eventWeight"), self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "(^GluGluHToTauTau.*125.*|^VBFHToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ggHEstimation(HTTEstimation):
    def __init__(self, era, directory, channel):
        super(HTTEstimation, self).__init__(
            name="ggH",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_files(self):
        query = {
            "process": "^GluGluHToTauTau.*125.*",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class qqHEstimation(HTTEstimation):
    def __init__(self, era, directory, channel):
        super(HTTEstimation, self).__init__(
            name="qqH",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_files(self):
        query = {
            "process": "^VBFHToTauTau.*125.*",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class VHEstimation(HTTEstimation):
    def __init__(self, era, directory, channel):
        super(HTTEstimation, self).__init__(
            name="VH",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_files(self):
        query = {
            "process": "(^W(minus|plus)HToTauTau.*125.*|^ZHToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ZTTEstimation(EstimationMethod):
    def __init__(self, era, directory, channel):
        super(ZTTEstimation, self).__init__(
            name="ZTT",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_weights(self):
        return Weights(
            Weight("eventWeight", "eventWeight"),
            Weight("zPtReweightWeight", "zPtReweightWeight"),
            Weight("((gen_match_2 == 5)*0.95 + (gen_match_2 != 5))",
                   "hadronic_tau_sf"),
            Weight(
                "((((genbosonmass >= 150.0 && (npartons == 0 || npartons >= 5))*3.95423374e-5) + ((genbosonmass >= 150.0 && npartons == 1)*1.27486147e-5) + ((genbosonmass >= 150.0 && npartons == 2)*1.3012785e-5) + ((genbosonmass >= 150.0 && npartons == 3)*1.33802133e-5) + ((genbosonmass >= 150.0 && npartons == 4)*1.09698723e-5)+((genbosonmass >= 50.0 && genbosonmass < 150.0 && (npartons == 0 || npartons >= 5))*3.95423374e-5) + ((genbosonmass >= 50.0 && genbosonmass < 150.0 && npartons == 1)*1.27486147e-5) + ((genbosonmass >= 50.0 && genbosonmass < 150.0 && npartons == 2)*1.3012785e-5) + ((genbosonmass >= 50.0 && genbosonmass < 150.0 && npartons == 3)*1.33802133e-5) + ((genbosonmass >= 50.0 && genbosonmass < 150.0 && npartons == 4)*1.09698723e-5)+((genbosonmass < 50.0)*numberGeneratedEventsWeight*crossSectionPerEventWeight))/(numberGeneratedEventsWeight*crossSectionPerEventWeight*sampleStitchingWeight))",
                "z_stitching_weight"), self.era.lumi_weight)

    def get_cuts(self):
        return Cuts(Cut("gen_match_2==5", "ztt_genmatch_mt"))

    def get_files(self):
        query = {
            "process":
            "(DYJetsToLL_M10to50|DYJetsToLL_M50|DY1JetsToLL_M50|DY2JetsToLL_M50|DY3JetsToLL_M50|DY4JetsToLL_M50)",
            "data":
            False,
            "campaign":
            self._mc_campaign,
            "generator":
            "madgraph\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ZLLEstimation(ZTTEstimation):
    def __init__(self, era, directory, channel):
        super(ZTTEstimation, self).__init__(
            name="ZLL",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("(gen_match_2<5||gen_match_2==6)", "zll_genmatch_mt"))

    def get_weights(self):
        ztt_weights = super(ZLLEstimation, self).get_weights()
        return ztt_weights + Weights(
            Weight(
                "(((decayMode_2 == 0)*1.0) + ((decayMode_2 == 1 || decayMode_2 == 2)*1.0) + ((decayMode_2 == 10)*1.0))",
                "decay_mode_reweight"))


class ZLLEstimationMTSM(ZLLEstimation):
    def get_weights(self):
        ztt_weights = super(ZLLEstimation, self).get_weights()
        return ztt_weights + Weights(
            Weight(
                "(((decayMode_2 == 0)*0.75) + ((decayMode_2 == 1 || decayMode_2 == 2)*1.0) + ((decayMode_2 == 10)*1.0))",
                "decay_mode_reweight"))


class ZLLEstimationETSM(ZLLEstimation):
    def get_weights(self):
        ztt_weights = super(ZLLEstimation, self).get_weights()
        return ztt_weights + Weights(
            Weight(
                "(((decayMode_2 == 0)*0.98) + ((decayMode_2 == 1 || decayMode_2 == 2)*1.2) + ((decayMode_2 == 10)*1.0))",
                "decay_mode_reweight"))


class ZLEstimationMT(ZTTEstimation):
    def __init__(self, era, directory, channel):
        super(ZTTEstimation, self).__init__(
            name="ZL",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("gen_match_2<5", "zl_genmatch_mt"))


class ZLEstimationMTSM(ZLEstimationMT):
    def get_weights(self):
        ztt_weights = super(ZLEstimationMT, self).get_weights()
        return ztt_weights + Weights(
            Weight(
                "(((decayMode_2 == 0)*0.75) + ((decayMode_2 == 1 || decayMode_2 == 2)*1.0) + ((decayMode_2 == 10)*1.0))",
                "decay_mode_reweight"))


class ZLEstimationETSM(ZLEstimationMT):
    def get_weights(self):
        ztt_weights = super(ZLEstimationMT, self).get_weights()
        return ztt_weights + Weights(
            Weight(
                "(((decayMode_2 == 0)*0.98) + ((decayMode_2 == 1 || decayMode_2 == 2)*1.2) + ((decayMode_2 == 10)*1.0))",
                "decay_mode_reweight"))


class ZJEstimationMT(ZTTEstimation):
    def __init__(self, era, directory, channel):
        super(ZTTEstimation, self).__init__(
            name="ZJ",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("gen_match_2==6", "zj_genmatch_mt"))


# et is equivalent to mt
class ZJEstimationET(ZJEstimationMT):
    pass


class ZLEstimationET(ZLEstimationMT):
    pass


# tt is equivalent to mt
class ZJEstimationTT(ZJEstimationMT):
    pass


class ZLEstimationTT(ZLEstimationMT):
    pass


class WEstimation(EstimationMethod):
    def __init__(self, era, directory, channel):
        super(WEstimation, self).__init__(
            name="W",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_weights(self):
        return Weights(
            Weight(
                "(((npartons == 0 || npartons >= 5)*7.09390278348407e-4) + ((npartons == 1)*1.90063898596475e-4) + ((npartons == 2)*5.8529964471165e-5) + ((npartons == 3)*1.9206444928444e-5) + ((npartons == 4)*1.923548021385e-5))/(numberGeneratedEventsWeight*crossSectionPerEventWeight*sampleStitchingWeight)",
                "wj_stitching_weight"),
            Weight("((gen_match_2 == 5)*0.95 + (gen_match_2 != 5))",
                   "hadronic_tau_sf"),
            Weight("eventWeight", "eventWeight"), self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "W.*JetsToLNu",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class TTEstimation(EstimationMethod):
    def __init__(self, era, directory, channel):
        super(TTEstimation, self).__init__(
            name="TT",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_weights(self):
        return Weights(
            Weight("topPtReweightWeight", "topPtReweightWeight"),
            Weight("eventWeight", "eventWeight"),
            Weight("((gen_match_2 == 5)*0.95 + (gen_match_2 != 5))",
                   "hadronic_tau_sf"), self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^TT$",
            "data": False,
            "campaign": self._mc_campaign
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class TTTEstimationMT(TTEstimation):
    def __init__(self, era, directory, channel):
        super(TTEstimation, self).__init__(
            name="TTT",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("gen_match_2==5", "ttt_genmatch_mt"))


class TTJEstimationMT(TTEstimation):
    def __init__(self, era, directory, channel):
        super(TTEstimation, self).__init__(
            name="TTJ",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("gen_match_2!=5", "ttj_genmatch_mt"))


class TTTEstimationET(TTTEstimationMT):
    pass


class TTJEstimationET(TTJEstimationMT):
    pass


class TTTEstimationTT(TTTEstimationMT):
    pass


class TTJEstimationTT(TTJEstimationMT):
    pass


class VVEstimation(EstimationMethod):
    def __init__(self, era, directory, channel):
        super(VVEstimation, self).__init__(
            name="VV",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_weights(self):
        return Weights(
            Weight("((gen_match_2 == 5)*0.95 + (gen_match_2 != 5))",
                   "hadronic_tau_sf"), Weight("eventWeight", "eventWeight"))

    def get_files(self):
        query = {
            "process":
            "(WWTo1L1Nu2Q|" + "WZTo1L1Nu2Q|" + "WZTo1L3Nu|" + "WZTo2L2Q|" +
            "ZZTo2L2Q" + ")",
            "data":
            False,
            "campaign":
            self._mc_campaign,
            "generator":
            "amcatnlo-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)

        query = {
            "process": "ZZTo4L",
            "extension": "ext1",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "amcatnlo-pythia8"
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)

        query = {
            "process": "WZJToLLLNu",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "pythia8"
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)

        query = {
            "process":
            "(STt-channelantitop4finclusiveDecays|STt-channeltop4finclusiveDecays|STtWantitop5finclusiveDecays|STtWtop5finclusiveDecays)",
            "data":
            False,
            "campaign":
            self._mc_campaign
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)

        log_query(self.name, "<optimzed out>", files)
        return self.artus_file_names(files)


class QCDEstimationET(SStoOSEstimationMethod):
    def __init__(self, era, directory, channel, bg_processes, data_process):
        super(QCDEstimationET, self).__init__(
            name="QCD",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            bg_processes=bg_processes,
            data_process=data_process)


class QCDEstimationMT(QCDEstimationET):
    pass


class QCDEstimationTT(ABCDEstimationMethod):
    def __init__(self, era, directory, channel, bg_processes, data_process):
        super(QCDEstimationTT, self).__init__(
            name="QCD",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            bg_processes=bg_processes,
            data_process=data_process,
            AC_cut_names=[ # cuts to be removed to include region for shape derivation
                "tau_2_iso"
            ],
            BD_cuts=[      # cuts to be applied to restrict to region for shape derivation
                Cut("byTightIsolationMVArun2v1DBoldDMwLT_2<0.5",
                    "tau_2_iso"),
                Cut("byLooseIsolationMVArun2v1DBoldDMwLT_2>0.5",
                    "tau_2_iso_loose")
            ],
            AB_cut_names=[ # cuts to be removed to include region for the determination of the extrapolation derivation
                "os"
            ],
            CD_cuts=[      # cuts to be applied to restrict to region for the determination of the extrapolation derivation
                Cut("q_1*q_2>0", "ss")
            ]
        )


class WEstimationWithQCD(EstimationMethod):
    def __init__(self, era, directory, channel, bg_processes, data_process, w_process, qcd_ss_to_os_extrapolation_factor):
        super(WEstimationWithQCD, self).__init__(
            name="WJets",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign=None)
        self._bg_processes = bg_processes
        self._data_process = data_process
        self._w_process = w_process
        self._qcd_ss_to_os_extrapolation_factor = qcd_ss_to_os_extrapolation_factor

    def create_root_objects(self, systematic):

        # create category for MC WJets shape estimation in the signal region
        signal_region = copy.deepcopy(systematic.category)
        signal_region.name = (signal_region.name+"_for_wjets_mc").lstrip(self.channel.name+"_")

        # create control regions for W yield estimation
        high_mt_ss_control_region = copy.deepcopy(systematic.category)
        high_mt_ss_control_region.name = "wjets_high_mt_ss_cr"
        high_mt_ss_control_region._variable = None

        high_mt_ss_control_region.cuts.remove("mt")
        high_mt_ss_control_region.cuts.remove("os")

        high_mt_ss_control_region.cuts.add(Cut("mt_1>70","mt"))
        high_mt_ss_control_region.cuts.add(Cut("q_1*q_2>0","ss"))

        # create control regions for W high mt to low mt extrapolation factor
        high_mt_os_control_region = copy.deepcopy(systematic.category) # this one also used for W yield estimation
        high_mt_os_control_region.name = "wjets_high_mt_os_cr"
        high_mt_os_control_region._variable = None

        high_mt_os_control_region.cuts.remove("mt")

        high_mt_os_control_region.cuts.add(Cut("mt_1>70","mt"))

        low_mt_os_control_region = copy.deepcopy(systematic.category)
        low_mt_os_control_region.name = "wjets_low_mt_os_cr"
        low_mt_os_control_region._variable = None

        low_mt_ss_control_region = copy.deepcopy(systematic.category)
        low_mt_ss_control_region.name = "wjets_low_mt_ss_cr"
        low_mt_ss_control_region._variable = None

        low_mt_ss_control_region.cuts.remove("os")

        low_mt_ss_control_region.cuts.add(Cut("q_1*q_2>0","ss"))

        # create control regions for W ss to os extrapolation factor
        inclusive_os_control_region = copy.deepcopy(systematic.category)
        inclusive_os_control_region.name = "wjets_os_cr"
        inclusive_os_control_region._variable = None

        inclusive_os_control_region.cuts.remove("mt")

        inclusive_ss_control_region = copy.deepcopy(systematic.category)
        inclusive_ss_control_region.name = "wjets_ss_cr"
        inclusive_ss_control_region._variable = None

        inclusive_ss_control_region.cuts.remove("mt")
        inclusive_ss_control_region.cuts.remove("os")

        inclusive_ss_control_region.cuts.add(Cut("q_1*q_2>0","ss"))

        # initialize root objects and systematics
        root_objects = []
        systematic._WandQCD_systematics = []

        # for extrapolation factors: only W MC is needed
        for category in [high_mt_os_control_region, low_mt_os_control_region, high_mt_ss_control_region, low_mt_ss_control_region, inclusive_os_control_region, inclusive_ss_control_region]:
            s = Systematic(
                category=category,
                process=self._w_process,
                analysis=systematic.analysis,
                era=self.era,
                variation=systematic.variation,
                mass=125)
            systematic._WandQCD_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects

        # for yields in high mt control regions: data and other bg processes needed
        for process in [self._data_process] + self._bg_processes:
            for category in [high_mt_os_control_region, high_mt_ss_control_region]:
                s = Systematic(
                    category=category,
                    process=process,
                    analysis=systematic.analysis,
                    era=self.era,
                    variation=systematic.variation,
                    mass=125)
                systematic._WandQCD_systematics.append(s)
                s.create_root_objects()
                root_objects += s.root_objects

        # for signal region shape
        s = Systematic(category=signal_region,process=self._w_process,analysis=systematic.analysis,era=self.era,variation=systematic.variation,mass=125)
        systematic._WandQCD_systematics.append(s)
        s.create_root_objects()
        root_objects += s.root_objects

        return root_objects

    def do_estimation(self, systematic):
        if not hasattr(systematic, "_WandQCD_systematics"):
            logger.fatal(
                "Systematic %s does not have attribute _WandQCD_systematics needed for WandQCD estimation.",
                systematic.name)
            raise Exception

        # Sort shapes and counts
        wjets_mc_shape = None
        wjets_high_mt_ss_cr_counts = {}
        wjets_high_mt_os_cr_counts = {}
        wjets_low_mt_os_cr_count = None
        wjets_low_mt_ss_cr_count = None
        wjets_os_cr_count = None
        wjets_ss_cr_count = None
        for s in systematic._WandQCD_systematics:
            s.do_estimation()
            if s.category.name.endswith("for_wjets_mc"):
                wjets_mc_shape = s.shape
            elif s.category.name.endswith("wjets_high_mt_ss_cr"):
                wjets_high_mt_ss_cr_counts[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_high_mt_os_cr"):
                wjets_high_mt_os_cr_counts[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_low_mt_os_cr"):
                wjets_low_mt_os_cr_count = s.shape
            elif s.category.name.endswith("wjets_low_mt_ss_cr"):
                wjets_low_mt_ss_cr_count = s.shape
            elif s.category.name.endswith("wjets_os_cr"):
                wjets_os_cr_count = s.shape
            elif s.category.name.endswith("wjets_ss_cr"):
                wjets_ss_cr_count = s.shape

        # Determine extrapolation factors
        R_ss_to_os = wjets_os_cr_count.result/wjets_ss_cr_count.result

        wjets_integral_low_mt_os = wjets_low_mt_os_cr_count.result
        wjets_integral_high_mt_os = wjets_high_mt_os_cr_counts.pop(self._w_process.name).result

        R_high_to_low_mt_os = wjets_integral_low_mt_os/wjets_integral_high_mt_os
        R_high_to_low_mt_ss = wjets_low_mt_ss_cr_count.result/wjets_high_mt_ss_cr_counts.pop(self._w_process.name).result
        print "SS to OS extrapolation factor:",R_ss_to_os
        print "high to low mt os extrapolation factor:",R_high_to_low_mt_os
        print "high to low mt ss extrapolation factor:",R_high_to_low_mt_ss

        # Determine yields in wjets CRs
        print "Data yield in ss high mt region:",wjets_high_mt_ss_cr_counts[self._data_process.name].result
        high_mt_ss_yield = wjets_high_mt_ss_cr_counts.pop(self._data_process.name).result - sum(
            [s.result for s in wjets_high_mt_ss_cr_counts.values()])
        sum_mc = sum([s.result for s in wjets_high_mt_ss_cr_counts.values()])
        print "MC yield to be subtracted:",sum_mc
        for name,s in wjets_high_mt_ss_cr_counts.items():
            print name,":",s.result/sum_mc
        print "yield in ss high mt region:",high_mt_ss_yield

        print "Data yield in os high mt region:",wjets_high_mt_os_cr_counts[self._data_process.name].result
        high_mt_os_yield = wjets_high_mt_os_cr_counts.pop(self._data_process.name).result - sum(
            [s.result for s in wjets_high_mt_os_cr_counts.values()])
        sum_mc = sum([s.result for s in wjets_high_mt_os_cr_counts.values()])
        print "MC yield to be subtracted:",sum_mc
        for name,s in wjets_high_mt_os_cr_counts.items():
            print name,":",s.result/sum_mc
        print "yield in os high mt region:",high_mt_os_yield

        # Derive and normalize final shape
        print "MC yield in signal region:",wjets_integral_low_mt_os
        sf = R_ss_to_os*(high_mt_os_yield - self._qcd_ss_to_os_extrapolation_factor*high_mt_ss_yield)/(R_ss_to_os-self._qcd_ss_to_os_extrapolation_factor)/wjets_integral_high_mt_os
        estimated_yield = R_high_to_low_mt_os*R_ss_to_os*(high_mt_os_yield - self._qcd_ss_to_os_extrapolation_factor*high_mt_ss_yield)/(R_ss_to_os-self._qcd_ss_to_os_extrapolation_factor)
        print "Estimated yield in signal region:",estimated_yield
        print "Scale wjets by",sf
        wjets_shape = copy.deepcopy(wjets_mc_shape)
        wjets_shape.result.Scale(sf)

        # Rename root object accordingly
        wjets_shape.name = systematic.name

        # Replace negative entries by zeros and renormalize shape
        wjets_shape.replace_negative_entries_and_renormalize(tolerance=0.05)

        return wjets_shape

    # Data-driven estimation, no associated files and weights
    def get_files(self):
        raise NotImplementedError

    def get_weights(self):
        raise NotImplementedError


class QCDEstimationWithW(EstimationMethod):
    def __init__(self, era, directory, channel, bg_processes, data_process, w_process, qcd_ss_to_os_extrapolation_factor):
        super(QCDEstimationWithW, self).__init__(
            name="QCD",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign=None)
        self._bg_processes = bg_processes
        self._data_process = data_process
        self._w_process = w_process
        self._qcd_ss_to_os_extrapolation_factor = qcd_ss_to_os_extrapolation_factor

    def create_root_objects(self, systematic):

        # create category for WJets and QCD shape estimation in the qcd control region
        qcd_control_region = copy.deepcopy(systematic.category)
        qcd_control_region.name = (qcd_control_region.name+"_ss_for_qcd").lstrip(self.channel.name+"_")

        qcd_control_region.cuts.remove("os")

        qcd_control_region.cuts.add(Cut("q_1*q_2>0","ss"))

        # create control regions for W yield estimation
        high_mt_ss_control_region = copy.deepcopy(systematic.category)
        high_mt_ss_control_region.name = "wjets_high_mt_ss_cr"
        high_mt_ss_control_region._variable = None

        high_mt_ss_control_region.cuts.remove("mt")
        high_mt_ss_control_region.cuts.remove("os")

        high_mt_ss_control_region.cuts.add(Cut("mt_1>70","mt"))
        high_mt_ss_control_region.cuts.add(Cut("q_1*q_2>0","ss"))

        # create control regions for W high mt to low mt extrapolation factor
        high_mt_os_control_region = copy.deepcopy(systematic.category) # this one also used for W yield estimation
        high_mt_os_control_region.name = "wjets_high_mt_os_cr"
        high_mt_os_control_region._variable = None

        high_mt_os_control_region.cuts.remove("mt")

        high_mt_os_control_region.cuts.add(Cut("mt_1>70","mt"))

        low_mt_os_control_region = copy.deepcopy(systematic.category)
        low_mt_os_control_region.name = "wjets_low_mt_os_cr"
        low_mt_os_control_region._variable = None

        low_mt_ss_control_region = copy.deepcopy(systematic.category)
        low_mt_ss_control_region.name = "wjets_low_mt_ss_cr"
        low_mt_ss_control_region._variable = None

        low_mt_ss_control_region.cuts.remove("os")

        low_mt_ss_control_region.cuts.add(Cut("q_1*q_2>0","ss"))

        # create control regions for W ss to os extrapolation factor
        inclusive_os_control_region = copy.deepcopy(systematic.category)
        inclusive_os_control_region.name = "wjets_os_cr"
        inclusive_os_control_region._variable = None

        inclusive_os_control_region.cuts.remove("mt")

        inclusive_ss_control_region = copy.deepcopy(systematic.category)
        inclusive_ss_control_region.name = "wjets_ss_cr"
        inclusive_ss_control_region._variable = None

        inclusive_ss_control_region.cuts.remove("mt")
        inclusive_ss_control_region.cuts.remove("os")

        inclusive_ss_control_region.cuts.add(Cut("q_1*q_2>0","ss"))

        # initialize root objects and systematics
        root_objects = []
        systematic._WandQCD_systematics = []

        # for extrapolation factors: only W MC is needed
        for category in [high_mt_os_control_region, low_mt_os_control_region, high_mt_ss_control_region, low_mt_ss_control_region, inclusive_os_control_region, inclusive_ss_control_region]:
            s = Systematic(
                category=category,
                process=self._w_process,
                analysis=systematic.analysis,
                era=self.era,
                variation=systematic.variation,
                mass=125)
            systematic._WandQCD_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects

        # for yields in high mt control regions: data and other bg processes needed
        for process in [self._data_process] + self._bg_processes:
            for category in [high_mt_os_control_region, high_mt_ss_control_region]:
                s = Systematic(
                    category=category,
                    process=process,
                    analysis=systematic.analysis,
                    era=self.era,
                    variation=systematic.variation,
                    mass=125)
                systematic._WandQCD_systematics.append(s)
                s.create_root_objects()
                root_objects += s.root_objects

        # for Wjets and QCD shape
        for process in [self._data_process,self._w_process] + self._bg_processes:
            s = Systematic(category=qcd_control_region,process=process,analysis=systematic.analysis,era=self.era,variation=systematic.variation,mass=125)
            systematic._WandQCD_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects

        return root_objects

    def do_estimation(self, systematic):
        if not hasattr(systematic, "_WandQCD_systematics"):
            logger.fatal(
                "Systematic %s does not have attribute _WandQCD_systematics needed for WandQCD estimation.",
                systematic.name)
            raise Exception

        # Sort shapes and counts
        qcd_control_region_shapes = {}
        wjets_high_mt_ss_cr_counts = {}
        wjets_high_mt_os_cr_counts = {}
        wjets_low_mt_os_cr_count = None
        wjets_low_mt_ss_cr_count = None
        wjets_os_cr_count = None
        wjets_ss_cr_count = None
        for s in systematic._WandQCD_systematics:
            s.do_estimation()
            if s.category.name.endswith("ss_for_qcd"):
                qcd_control_region_shapes[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_high_mt_ss_cr"):
                wjets_high_mt_ss_cr_counts[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_high_mt_os_cr"):
                wjets_high_mt_os_cr_counts[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_low_mt_os_cr"):
                wjets_low_mt_os_cr_count = s.shape
            elif s.category.name.endswith("wjets_low_mt_ss_cr"):
                wjets_low_mt_ss_cr_count = s.shape
            elif s.category.name.endswith("wjets_os_cr"):
                wjets_os_cr_count = s.shape
            elif s.category.name.endswith("wjets_ss_cr"):
                wjets_ss_cr_count = s.shape

        # Determine extrapolation factors
        R_ss_to_os = wjets_os_cr_count.result/wjets_ss_cr_count.result

        wjets_integral_low_mt_ss = wjets_low_mt_ss_cr_count.result
        wjets_integral_high_mt_ss = wjets_high_mt_ss_cr_counts.pop(self._w_process.name).result

        R_high_to_low_mt_os = wjets_low_mt_os_cr_count.result/wjets_high_mt_os_cr_counts.pop(self._w_process.name).result
        R_high_to_low_mt_ss = wjets_integral_low_mt_ss/wjets_integral_high_mt_ss
        print "SS to OS extrapolation factor:",R_ss_to_os
        print "high to low mt os extrapolation factor:",R_high_to_low_mt_os
        print "high to low mt ss extrapolation factor:",R_high_to_low_mt_ss

        # Determine yields in wjets CRs
        high_mt_ss_yield = wjets_high_mt_ss_cr_counts.pop(self._data_process.name).result - sum(
            [s.result for s in wjets_high_mt_ss_cr_counts.values()])

        high_mt_os_yield = wjets_high_mt_os_cr_counts.pop(self._data_process.name).result - sum(
            [s.result for s in wjets_high_mt_os_cr_counts.values()])

        print "yield in ss high mt region:",high_mt_ss_yield
        print "yield in os high mt region:",high_mt_os_yield

        # Derive and normalize final shape for QCD
        wjets_shape = qcd_control_region_shapes.pop(self._w_process.name)
        print "MC yield in qcd control region for wjets:",wjets_integral_low_mt_ss
        sf = (high_mt_os_yield - self._qcd_ss_to_os_extrapolation_factor*high_mt_ss_yield)/(R_ss_to_os-self._qcd_ss_to_os_extrapolation_factor)/wjets_integral_high_mt_ss
        estimated_yield = R_high_to_low_mt_ss*(high_mt_os_yield - self._qcd_ss_to_os_extrapolation_factor*high_mt_ss_yield)/(R_ss_to_os-self._qcd_ss_to_os_extrapolation_factor)
        print "Estimated yield in qcd control region for wjets:",estimated_yield
        print "Scale wjets by",sf
        wjets_shape.result.Scale(sf)

        qcd_shape = copy.deepcopy(qcd_control_region_shapes.pop(self._data_process.name))
        qcd_shape.result.Add(wjets_shape.result,-1.0)
        for sh in qcd_control_region_shapes.values():
            qcd_shape.result.Add(sh.result,-1.0)
        qcd_shape.result.Scale(self._qcd_ss_to_os_extrapolation_factor)

        # Rename root object accordingly
        qcd_shape.name = systematic.name

        # Replace negative entries by zeros and renormalize shape
        qcd_shape.replace_negative_entries_and_renormalize(tolerance=0.05)

        return qcd_shape

    # Data-driven estimation, no associated files and weights
    def get_files(self):
        raise NotImplementedError

    def get_weights(self):
        raise NotImplementedError
