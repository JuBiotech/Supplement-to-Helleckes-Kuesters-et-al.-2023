import logging
import typing

import calibr8
import numpy
import pandas
import pymc
import aesara.tensor as at


_log = logging.getLogger(__file__)


class NADHFluorescenceModel(calibr8.BaseExponentialModelN):
    def __init__(self):
        super().__init__(
            independent_key="NADH_concentration",
            dependent_key="fluorescence",
            sigma_degree=1
        )


def _add_or_assert_coords(
    coords: typing.Dict[str, typing.Sequence], pmodel: pymc.Model
):
    """Ensures that the coords are available in the model."""
    for cname, cvalues in coords.items():
        if cname in pmodel.coords:
            numpy.testing.assert_array_equal(pmodel.coords[cname], cvalues)
        else:
            pmodel.add_coord(name=cname, values=cvalues)

class CatIBModel:
    def __init__(
        self,
        *,
        df_kinetics: pandas.DataFrame,
        cm_nadh: typing.Optional[calibr8.CalibrationModel],
        sd_assay: float = 0.1,
    ):
        """Create a model for the analysis of NADH fluorescence of GDH-CatIBs in enzymatic assay.

        Parameters
        ----------
        df_kinetics : pandas.DataFrame
            Observations of the GDH assay.
            Index: [strain, kinetic_id, concentration_factor, cycle]
            Columns: [time, value]
        cm_nadh : calibr8.CalibrationModel
            The calibration model for 4-nitrophenol absorbance/concentration
        sd_assay : float
            Approximate pipetting error from input well to fluorescence assay well.
            0.1 ≙ 10 % relative error
        """
        self.df_kinetics = df_kinetics
        self.cm_nadh = cm_nadh
        self.idx_KIDtoS = []
        self.idx_KIDtoCID = []
        self.idx_KIDtoRun = []
        self.sd_assay=sd_assay
        self.strains = sorted([
            'RoSa126', 'RoSa124','RoSa123','RoSa122', 'RoSa12', 'RoSa11', 'RoSa10', 'RoSa9',
            'RoSa121', 'RoSa108', 'RoSa107', 'RoSa106','RoSa32', 'RoSa31', 'RoSa30', 'RoSa29',
            'KK166', 'RoSa225', 'RoSa224', 'RoSa223', 'RoSa105', 'RoSa104', 'RoSa103', 'RoSa102',
            'RoSa28', 'RoSa27', 'RoSa26', 'RoSa25', 'RoSa221', 'RoSa220', 'RoSa219', 'RoSa217',
            'RoSa101', 'RoSa100', 'RoSa99', 'RoSa98', 'KK184', 'KK183', 'KK181', 'KK179',
            'RoSa216',  'RoSa215', 'RoSa214', 'RoSa132', 'RoSa97', 'RoSa96', 'RoSa95', 'RoSa94',
            'KK177', 'KK176', 'KK174', 'KK173', 'RoSa131', 'RoSa130', 'RoSa128', 'RoSa127',
            'RoSa52', 'RoSa16', 'RoSa15', 'RoSa13', 'KK172', 'KK171', 'KK169','KK167'
        ])
        self.kinetic_ids = None
        self.column_ids = None

        with pymc.Model() as self.pmodel:
            self._model_concentrations()
            self._model_nadh()
        super().__init__()

    def _model_concentrations(self):

        # first check that the input dataframe corresponds to the standard format:

        # input checks successful! Can now start building the model...
        # First get a handle on the pymc model
        pmodel = pymc.modelcontext(None)

        S = "strain"
        KID = "kinetic_id"
        CID = "column_id"
        R = "run"
        
        self.kinetic_ids = list(sorted(set(self.df_kinetics.index.get_level_values(KID))))
        self.column_ids = list(sorted(set(self.df_kinetics.index.get_level_values(CID))))
        self.runs = list(sorted(set(self.df_kinetics.index.get_level_values(R))))

        # The model describes the data in the "long form". It needs some lists for indexing to map
        # between variables of different dimensionality.
        # For example the `idx_KIDtoS` maps each kinetic id (KID) to the corresponding strain (S).
        concentration_factors = []
        for kid in self.kinetic_ids:
            # make a mapping to generate data in the dim of KID from the list of strains
            strain = numpy.unique(self.df_kinetics.loc[:,kid,:,:,:].index.get_level_values(S))
            assert len(strain==1), "Kinetic id should be unique and only belong to one strain."
            self.idx_KIDtoS.append(self.strains.index(strain[0]))

            # make a mapping to generate data in the dim of KID from the list of columns
            column_id = numpy.unique(self.df_kinetics.loc[:,kid,:,:,:].index.get_level_values(CID))
            assert len(column_id==1), "Kinetic id should be unique and only belong to one column."
            self.idx_KIDtoCID.append(self.column_ids.index(column_id[0]))

            # make a mapping to generate data in the dim of KID from the run_ids
            run = numpy.unique(self.df_kinetics.loc[:,kid,:,:,:].index.get_level_values(R))
            assert len(run==1), "Kinetic id should be unique and only belong to one run."
            self.idx_KIDtoRun.append(self.runs.index(run[0]))

            # we also extract the concentration factor for each kinetic id for a later part of the model
            cf = numpy.unique(self.df_kinetics.loc[:,kid,:,:,:].index.get_level_values("concentration_factor").to_numpy(dtype=float))
            assert len(cf==1), "The concentration factor should be the same for all cycles of one kinetic_id."
            concentration_factors.append(cf[0])
        assert numpy.shape(self.idx_KIDtoS) == (len(self.kinetic_ids),), \
        f"{numpy.shape(self.idx_KIDtoS)},{(len(self.kinetic_ids),)}"
        assert max(self.idx_KIDtoS) == len(self.strains) - 1
        assert numpy.shape(self.idx_KIDtoCID) == (len(self.kinetic_ids),), \
        f"{numpy.shape(self.idx_KIDtoCID)},{(len(self.kinetic_ids),)}"
        assert max(self.idx_KIDtoCID) == 7
        assert numpy.shape(self.idx_KIDtoRun) == (len(self.kinetic_ids),), \
        f"{numpy.shape(self.idx_KIDtoRun)},{(len(self.kinetic_ids),)}"
        assert max(self.idx_KIDtoRun) == len(self.runs) - 1
        # Now create "coords" that describe the relevant dimensions
        _add_or_assert_coords(            
            {
                S: numpy.array(self.strains, dtype=str),
                KID: numpy.array(self.kinetic_ids, dtype=str),
                CID: numpy.array(self.column_ids, dtype=int),
                R: numpy.array(self.runs, dtype=str)
            },
            pmodel,
        )
        cf_input = pymc.Data(
            "cf_input", numpy.array(concentration_factors), dims=(KID,)
        )
        pymc.Lognormal(
            "cf_nadh_assay",
            mu=at.log(cf_input),
            sd=self.sd_assay,
            dims=(KID,),
        )

        return

    def _model_nadh(self):
        t_obs = numpy.vstack([
            self.df_kinetics.xs(kid, level=('kinetic_id')).time 
            for kid in self.kinetic_ids
        ])
        y_obs = numpy.vstack([
            self.df_kinetics.xs(kid, level=('kinetic_id')).value 
            for kid in self.kinetic_ids
        ])

        pmodel = pymc.modelcontext(None)
        _add_or_assert_coords(
            {
                "assay_cycle": numpy.arange(t_obs.shape[1]),
            },
            pmodel,
        )
        cf_nadh_assay = pmodel["cf_nadh_assay"]
        # Again, we assume one undiluted activity and calculate the activities in the sample/replicate wells
        # deterministically from the relative concentrations above.
        k_sd= pymc.HalfNormal("k_sd", sd=0.2) 
        k_mean = pymc.HalfNormal("k_mean", sd=0.1)
        k = pymc.Lognormal("k", mu=at.log(k_mean), sd=k_sd, dims=("strain",), )
        batch_effect = pymc.Lognormal(
            "batch_effect",
            mu=0,
            sd=0.1,
            dims=("run")
        )
        k_batch = pymc.Deterministic(
            "k_batch",
            k[self.idx_KIDtoS] * batch_effect[self.idx_KIDtoRun],
            dims=("kinetic_id")
        )
        # The inputs are diluted 50x into the assay.
        dilution_factor = pymc.Data("dilution_factor", 50)
        k_assay = pymc.Deterministic(
            "k_assay",
            cf_nadh_assay * k_batch / dilution_factor, #cf_nadh_assay brings the uncertainty of the dilution
            dims=("kinetic_id",),
        )

        # prediction
        cutinase_time = pymc.Data(
            "assay_time", t_obs, dims=("kinetic_id", "assay_cycle")
        )
        t_offset = pymc.Lognormal(
            "t_offset",
            mu=at.log(0.283),
            sd=0.01,
        )
        t_plate_to_reader = pymc.Lognormal(
            "t_plate_to_reader",
            mu=at.log(0.5),
            sd=0.05,
        )
        t_column = pymc.Deterministic(
            "t_column",
            numpy.arange(len(self.column_ids)) * t_offset + t_plate_to_reader,
            dims=("column_id")
        )
        S0 = pymc.Lognormal(
            "S0", 
            mu=at.log(0.2), 
            sd=0.1
        ) #truth should be 0.4 mM per Well (250 uL im Well)

        product_concentration = pymc.Deterministic(
            "product_concentration",
            S0 * (1 - at.exp(-k_assay[:, None] * (cutinase_time + t_column[self.idx_KIDtoCID, None]))),
            dims=("kinetic_id", "assay_cycle"),
        )

        # The reaction product (NADH) is measured in [mmol/L] = [µmol/mL] via the error model.
        # Our prediction is also in µmol/mL.
        fluorescence = pymc.Data(
            "fluorescence", y_obs , dims=("kinetic_id", "assay_cycle")
        )

        log_likelihood = self.cm_nadh.loglikelihood(
            y=fluorescence,
            x=product_concentration,
            name="nadh_all",
            dims=("kinetic_id", "assay_cycle"),
            dependent_key=self.cm_nadh.dependent_key,
        )
        return

 
    def summary(self):
        for cname, cvals in self.pmodel.coords.items():
            shape = numpy.shape(cvals)
            vals = numpy.asarray(cvals).flatten()
            if numpy.prod(shape) < 5:
                examples = ", ".join([str(v) for v in vals])
            else:
                examples = ", ".join([str(v) for v in vals][:3])
                examples += f", …, {cvals[-1]}"
            print(f"{cname: <20}{shape}\t{examples}")
