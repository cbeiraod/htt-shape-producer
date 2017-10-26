# -*- coding: utf-8 -*-

import copy
import os

from estimation_methods import EstimationMethod
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


# TODO: Does this work?
class HttEstimation(EstimationMethod):
    def __init__(self, era, directory, channel):
        super(HttEstimation, self).__init__(
            name="Htt",
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


class ZttEstimation(EstimationMethod):
    def __init__(self, era, directory, channel):
        super(ZttEstimation, self).__init__(
            name="Ztt",
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
        return Cuts(Cut("gen_match_2==5",
                        "ztt_genmatch_mt"))  # FIXME: Doubles with weights?

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


class ZllEstimation(ZttEstimation):
    def __init__(self, era, directory, channel):
        super(ZttEstimation, self).__init__(
            name="Zll",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("(gen_match_2<5||gen_match_2==6)", "zll_genmatch_mt"))

    def get_weights(self):
        ztt_weights = super(ZllEstimation, self).get_weights()
        return ztt_weights + Weights(
            Weight(
                "(((decayMode_2 == 0)*1.0) + ((decayMode_2 == 1 || decayMode_2 == 2)*1.0) + ((decayMode_2 == 10)*1.0))",
                "decay_mode_reweight"))


class WJetsEstimation(EstimationMethod):
    def __init__(self, era, directory, channel):
        super(WJetsEstimation, self).__init__(
            name="WJets",
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
            name="tt",
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


class QCDEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, bg_processes, data_process):
        super(QCDEstimation, self).__init__(
            name="QCD",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign=None)
        self._bg_processes = [copy.deepcopy(p) for p in bg_processes]
        self._data_process = copy.deepcopy(data_process)

    def create_root_objects(self, systematic):
        ss_category = copy.deepcopy(systematic.category)
        ss_category.cuts.get("os").name = "ss"
        ss_category.cuts.get("ss").invert()
        ss_category.name = ss_category._name + "_ss"

        root_objects = []
        systematic._qcd_systematics = []
        for process in [self._data_process] + self._bg_processes:
            s = Systematic(
                category=ss_category,
                process=process,
                analysis=systematic.analysis,
                era=self.era,
                variation=systematic.variation,
                mass=125)
            systematic._qcd_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects
        return root_objects

    def do_estimation(self, systematic):
        if not hasattr(systematic, "_qcd_systematics"):
            logger.fatal(
                "Systematic %s does not have attribute _qcd_systematics needed for QCD estimation.",
                systematic.name)
            raise Exception

        # Create shapes
        for s in systematic._qcd_systematics:
            s.do_estimation()

        # Data shape
        shape = systematic._qcd_systematics[0].shape

        # Subtract MC shapes from data shape
        for s in systematic._qcd_systematics[1:]:
            shape.result.Add(s.shape.result, -1.0)

        # Test that not a single bin in TH1F shape.result is negative
        if shape.has_negative_entries():
            logger.fatal(
                "Subtraction of Monte Carlo from data results in negative number of events."
            )
            raise Exception

        # Rename root object accordingly
        shape.name = systematic.name
        return shape

    # Data-driven estimation, no associated files and weights
    def get_files(self):
        raise NotImplementedError

    def get_weights(self):
        raise NotImplementedError



class QCDEstimationABCD(EstimationMethod):
    def __init__(self, era, directory, channel, bg_processes, data_process):
        super(QCDEstimationABCD, self).__init__(
            name="QCD",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign=None)
        self._bg_processes = [copy.deepcopy(p) for p in bg_processes]
        self._data_process = copy.deepcopy(data_process)

    def create_root_objects(self, systematic):
        ss_isolated_category = copy.deepcopy(systematic.category)
        ss_isolated_category.cuts.get("os").name = "ss"
        ss_isolated_category.cuts.get("ss").invert()
        ss_isolated_category.name = ss_isolated_category._name + "_ss_isolated"
        ss_isolated_category.variable._name += "_count"

        ss_antiisolated_category = copy.deepcopy(systematic.category)
        ss_antiisolated_category.cuts.get("os").name = "ss"
        ss_antiisolated_category.cuts.get("ss").invert()
        ss_antiisolated_category.cuts.get("tau_2_iso").invert()
        ss_antiisolated_category.cuts.add(Cut("byLooseIsolationMVArun2v1DBoldDMwLT_2>0.5","tau_2_iso_loose"))
        ss_antiisolated_category.name = ss_antiisolated_category._name + "_ss_antiisolated"
        ss_antiisolated_category.variable._name += "_count"

        os_antiisolated_category = copy.deepcopy(systematic.category)
        os_antiisolated_category.cuts.get("tau_2_iso").invert()
        os_antiisolated_category.cuts.add(Cut("byLooseIsolationMVArun2v1DBoldDMwLT_2>0.5","tau_2_iso_loose"))
        os_antiisolated_category.name = os_antiisolated_category._name + "_os_antiisolated"

        root_objects = []
        systematic._qcd_systematics = []
        for process in [self._data_process] + self._bg_processes:
            for category in [ss_isolated_category, ss_antiisolated_category, os_antiisolated_category]:
                s = Systematic(
                    category=category,
                    process=process,
                    analysis=systematic.analysis,
                    era=self.era,
                    variation=systematic.variation,
                    mass=125)
                systematic._qcd_systematics.append(s)
                s.create_root_objects()
                root_objects += s.root_objects
        return root_objects

    def do_estimation(self, systematic):
        if not hasattr(systematic, "_qcd_systematics"):
            logger.fatal(
                "Systematic %s does not have attribute _qcd_systematics needed for QCD estimation.",
                systematic.name)
            raise Exception

        # Create shapes
        ss_isolated_shapes = {}
        ss_antiisolated_shapes = {}
        os_antiisolated_shapes = {}
        for s in systematic._qcd_systematics:
            s.do_estimation()
            if "ss_isolated" in s.category.name:
                ss_isolated_shapes[s.process.name] = s.shape
            elif "ss_antiisolated" in s.category.name:
                ss_antiisolated_shapes[s.process.name] = s.shape
            elif "os_antiisolated" in s.category.name:
                os_antiisolated_shapes[s.process.name] = s.shape

        print "amount of shapes: ss iso =",len(ss_isolated_shapes),"ss anti-iso =",len(ss_antiisolated_shapes),"os anti-iso =",len(os_antiisolated_shapes)
        for s in ss_isolated_shapes.values():
            print s.name,":",type(s)
        for s in ss_antiisolated_shapes.values():
            print s.name,":",type(s)
        for s in os_antiisolated_shapes.values():
            print s.name,":",type(s)

        ss_isolated_yield = ss_isolated_shapes.pop("data_obs").result - sum([s.result for s in ss_isolated_shapes.values()])
        ss_antiisolated_yield = ss_antiisolated_shapes.pop("data_obs").result - sum([s.result for s in ss_antiisolated_shapes.values()])

        iso_extrapolation_factor = ss_isolated_yield/ss_antiisolated_yield

        print "extrapolation factor from antiisolated to the isolated region:",iso_extrapolation_factor

        # Subtract MC shapes from data shape
        for s in os_antiisolated_shapes:
            if s != "data_obs":
                os_antiisolated_shapes["data_obs"].result.Add(os_antiisolated_shapes[s].result, -1.0)

        os_antiisolated_shapes["data_obs"].result.Scale(iso_extrapolation_factor)

        # Test that not a single bin in TH1F shape.result is negative
        if os_antiisolated_shapes["data_obs"].has_negative_entries():
            logger.fatal(
                "Subtraction of Monte Carlo from data results in negative number of events."
            )
            raise Exception

        # Rename root object accordingly
        os_antiisolated_shapes["data_obs"].name = systematic.name
        return os_antiisolated_shapes["data_obs"]

    # Data-driven estimation, no associated files and weights
    def get_files(self):
        raise NotImplementedError

    def get_weights(self):
        raise NotImplementedError
