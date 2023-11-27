import logging
import os
import pathlib
import typing

import numpy
import pandas


_log = logging.getLogger(__file__)


VARIANT_MAPPER = {
    "KK166": ("$Bs$"+"GDH-SG-L12KD"),
    "KK167": ("$Bs$"+ "GDH-SG-L24KD"),
    "KK169": ("$Bs$"+ "GDH-SG-L48KD"),
    "KK171": ("$Bs$"+ "GDH-SG-" + "$\mathrm{(L6KD)}_2$"),
    "KK172": ("$Bs$"+ "GDH-SG-" + "$\mathrm{(L6KD)}_4$"),
    "KK173": ("$Bs$"+ "GDH-SG-" + "$\mathrm{(L6KD)}_6$"),
    "KK174": ("$Bs$"+ "GDH-SG-" + "$\mathrm{(L6KD)}_8$"),
    "KK176": ("L12KD-SG-$Bs$GDH"),
    "KK177": ("L24KD-SG-$Bs$GDH"),
    "KK179": ("L48KD-SG-$Bs$GDH"),
    "KK181": ("$\mathrm{(L6KD)}_2$" + "-SG-$Bs$GDH"),
    "KK183": ("$\mathrm{(L6KD)}_6$" + "-SG-$Bs$GDH"),
    "KK184": ("$\mathrm{(L6KD)}_8$" + "-SG-$Bs$GDH"),
    "RoSa9": ("TDoT-PT-$Bs$GDH"),
    "RoSa10": ("TorA-PT-$Bs$GDH"),
    "RoSa11": ("ELK16-PT-$Bs$GDH"),
    "RoSa12": ("18AWT-PT-$Bs$GDH"),
    "RoSa13": ("L6KD-PT-$Bs$GDH"),
    "RoSa15": ("3HAMP-PT-$Bs$GDH"),
    "RoSa16": ("CBDCell-PT-$Bs$GDH"),
    "RoSa25": ("TDoT-SG-$Bs$GDH"),
    "RoSa26": ("TorA-SG-$Bs$GDH"),
    "RoSa27": ("ELK16-SG-$Bs$GDH"),
    "RoSa28": ("18AWT-SG-$Bs$GDH"),
    "RoSa29": ("L6KD-SG-$Bs$GDH"),
    "RoSa30": ("GFIL8-SG-$Bs$GDH"),
    "RoSa31": ("3HAMP-SG-$Bs$GDH"),
    "RoSa32": ("CBDCell-SG-$Bs$GDH"),
    "RoSa52": ("soluble $Bs$GDH"),
    "RoSa94": ("$Bs$"+ "GDH-PT-TorA"),
    "RoSa95": ("$Bs$"+ "GDH-PT-CBDCell"),
    "RoSa96": ("$Bs$"+ "GDH-PT-3HAMP"),
    "RoSa97": ("$Bs$"+ "GDH-PT-GFIL8"),
    "RoSa98": ("$Bs$"+ "GDH-PT-L6KD"),
    "RoSa99": ("$Bs$"+ "GDH-PT-18AWT"),
    "RoSa100": ("$Bs$"+ "GDH-PT-TDoT"),
    "RoSa101": ("$Bs$"+ "GDH-SG-ELK16"),
    "RoSa102": ("$Bs$"+ "GDH-SG-TorA"),
    "RoSa103": ("$Bs$"+ "GDH-SG-CBDCell"),
    "RoSa104": ("$Bs$"+ "GDH-SG-3HAMP"),
    "RoSa105": ("$Bs$"+ "GDH-SG-GFIL8"),
    "RoSa106": ("$Bs$"+ "GDH-SG-L6KD"),
    "RoSa107": ("$Bs$"+ "GDH-SG-18AWT"),
    "RoSa108": ("$Bs$"+ "GDH-SG-TDoT"),
    "RoSa121": ("$Bs$"+ "GDH-G1-L6KD"),
    "RoSa122": ("$Bs$"+ "GDH-P1-L6KD"),
    "RoSa123": ("$Bs$"+ "GDH-G2-L6KD"),
    "RoSa124": ("$Bs$"+ "GDH-P2-L6KD"),
    "RoSa126": ("$Bs$"+ "GDH-P3-L6KD"),
    "RoSa127": ("$Bs$"+ "GDH-G4-L6KD"),
    "RoSa128": ("$Bs$"+ "GDH-P4-L6KD"),
    "RoSa130": ("$Bs$"+ "GDH-P5-L6KD"),
    "RoSa131": ("$Bs$"+ "GDH-G10-L6KD"),
    "RoSa132": ("$Bs$"+ "GDH-P10-L6KD"),
    "RoSa214": ("L6KD-G1-$Bs$GDH"),
    "RoSa215": ("L6KD-P1-$Bs$GDH"),
    "RoSa216": ("L6KD-G2-$Bs$GDH"),
    "RoSa217": ("L6KD-P2-$Bs$GDH"),
    "RoSa219": ("L6KD-P3-$Bs$GDH"),
    "RoSa220": ("L6KD-G4-$Bs$GDH"),
    "RoSa221": ("L6KD-P4-$Bs$GDH"),
    "RoSa223": ("L6KD-P5-$Bs$GDH"),
    "RoSa224": ("L6KD-G10-$Bs$GDH"),
    "RoSa225": ("L6KD-P10-$Bs$GDH"),
}


def extract_fps_for_round(runinfo_df, fp_experiment, round_no=0):
    """Extract filepaths for pre-, main culture and assay for given round.
    
    Parameters
    ----------
    fp_experiment : os.PathLike
        Path to the experiment folder
    runinfo_df : pandas.DataFrame
        DataFrame containing run_ids and result file names for each round.
    round_no : int
        Number of round to parse files (Starting from 0)
        
    Returns
    -------
    fp_VK : pathlib.Path
        Filepath to preculture csv   
    fp_HK : pathlib.Path
        Filepath to main culture csv   
    fp_assay : pathlib.Path
        Filepath to fluorescence xml
    fp_layout : pathlib.Path
        Filepath to file 'screening{round_no}.xlsx'
        
    """
    assert runinfo_df.index.name == "sround", \
    f"Index column must be names 'sround', but is {runinfo_df.index.name}." 
    
    sub_df = runinfo_df[runinfo_df.index == round_no]
    fp_VK = fp_experiment / "Screening_1-Preculture" / \
            sub_df["run_id-VK"][round_no] / sub_df["result_file-VK"][round_no]
    fp_HK = fp_experiment / "Screening_2-Mainculture" / \
            sub_df["run_id-HK"][round_no] / sub_df["result_file-HK"][round_no]
    fp_assay = fp_experiment / "Screening_3-Purification_and_Assay" / \
               sub_df["run_id-Assay"][round_no] / sub_df["result_file(s)-Assay"][round_no].split(",")[0]
    fp_layout = fp_experiment / "Screening_1-Preculture" / \
            sub_df["run_id-VK"][round_no] / f"screening{round_no+1}.xlsx"
    if not fp_VK.exists():
        _log.warning(
            f"The preculture file was not found in expected location {fp_VK}.\
            Check file location and runinfo_df for correct information."
        )
    if not fp_HK.exists():
        _log.warning(
            f"The main culture file was not found in expected location {fp_HK}. \
            Check file location and runinfo_df for correct information."
        )
    if not fp_assay.exists():
        _log.warning(
            f"The assay data file was not found in expected location {fp_assay}. \
            Check file location and runinfo_df for correct information."
        )
    if not fp_layout.exists():
        _log.warning(
            f"A file 'screening{round_no}.xlsx' containing the information on sample and \
            standards layout should be placed in the preculture folder. Check the \
            file location."
        )
    
    return fp_VK, fp_HK, fp_assay, fp_layout


def read_fluorescence(
    filepath : os.PathLike,
    wells : typing.Optional[numpy.ndarray] = None
) -> pandas.DataFrame:
    """Parse NADH fluorescence data into DataFrame.
    
    Parameters
    ----------
    filepath : os.PathLike
        The path to the xml file with fluorescence readouts.
    wells : numpr.ndarray
        The assay wells to filter for.

    Returns
    -------
    df_fluorescence : pandas.DataFrame
        DataFrame with time and readout values.
    """
    df_fluorescence = pandas.read_excel(str(filepath).replace("xml", "xlsx"), index_col=[0,1])
    if wells is not None:
        df_fluorescence = df_fluorescence[df_fluorescence.index.get_level_values('well').isin(wells)]
    return df_fluorescence


def read_NADH_standards(fp_layout:pathlib.Path, fp_assay:pathlib.Path) -> pandas.DataFrame:
    """Parse NADH calibration data into DataFrame for calibr8 model.
    
    Parameters
    ----------
    fp_layout : pathlib.Path
        File path to the layout excel sheet ('screening1.xlsx').
    fp_assay : pathlib.Path
        File path to the assay xml data.
    
    Returns
    -------
    df_calibration: pandas.DataFrame
        DataFrame with concentrations and fluorescence values of NADH standards
    """
    
    layout_df = pandas.read_excel(fp_layout, sheet_name="standards")
    std_wells = layout_df.mtp_well.values
    std_concentrations = layout_df.concentration.values
    standards = read_fluorescence(filepath=fp_assay, wells=std_wells)
    # only use measurements from 30 cycle on and reduce replicates
    cycles = numpy.unique(standards.index.get_level_values('cycle'))[30::10]
    # concentrations of replicates, ordered Fortran style:
    readouts = standards[standards.index.get_level_values('cycle').isin(cycles)].value.values
    # numpy magic to expand it into the same shape as the wells:
    std_concentrations = numpy.repeat(
        std_concentrations, repeats=len(cycles)
    ) #dilution factor already included in excel sheet

    # compile into a nice DataFrame
    df_calibration = pandas.DataFrame(
        columns=["concentration", "fluorescence"]
    )
    df_calibration["concentration"] = std_concentrations
    df_calibration["fluorescence"] = readouts
    df_calibration.sort_values("concentration", inplace=True)
    return df_calibration


def read_kinetics_run(*, runinfo_df:pandas.DataFrame, fp_experiment:pathlib.Path, round_no:int=0, concentration_factor:int=1):
    """Reads the NADH fluorescence data for one run into a DataFrame.
    
    Parameters
    ----------
    runinfo_df : pandas.Dataframe
        index is the number of rounda, necessary columns are run_id-VK, result_file-VK, run_id-HK, result_file-HK, run_id-Assay, result_file(s)-Assay
    fp_experiment : pathlib.Path
        The path to the top-level folder containing all run folders.
    round_no : int
        current screening round (>=0, ascending number of rounds)
    concentration_factor : int
        default 1, can be adapted if different between runs

    Returns
    -------
    df_kinetics : pandas.DataFrame
        Index: 
            strain : unique ID for each strain
            kinetic_id : ID of format well_run
            run : Hagelkorn ID of the assay run
            column_id : 1-12, corresponding to MTP column
            concentration_factor : default 1, can be adapted if different between runs
            cycle : measurement cycle in plate reader
        Columns:
            time : time after start of Reader measurement (attention if pipetting of columns has a time shift)
            value : fluorescence readout
    """
    _, _, fp_assay, fp_layout = extract_fps_for_round(runinfo_df, fp_experiment, round_no)
    assay_run_id = fp_assay.parts[-2]
    # read the layout
    layout_df = pandas.read_excel(fp_layout, sheet_name="samples")
    sample_wells = layout_df.mtp_well.values
    strains = layout_df.strain.values
    strain_mapping = {well : strain for well, strain in zip(sample_wells, strains)}
    column_mapping = {well : well.split("0")[-1] for well in sample_wells}
    
    # read the kinetic data
    df_kinetics = read_fluorescence(fp_assay, wells=sample_wells)
    df_kinetics["strain"] = df_kinetics.index.get_level_values("well").map(strain_mapping)
    df_kinetics["column_id"] = df_kinetics.index.get_level_values("well").map(column_mapping)
    # add run name to the well for unique identifiers
    df_kinetics.rename(index={
        old_name: f"{old_name}_{assay_run_id}"
        for old_name in df_kinetics.index.get_level_values("well")
    }, inplace=True)
    df_kinetics.reset_index(inplace=True)
    df_kinetics["concentration_factor"] = concentration_factor
    df_kinetics["run"] = assay_run_id
    df_kinetics.set_index(["strain", "well", "run", "column_id", "concentration_factor", "cycle"], inplace=True)
    df_kinetics.index.rename(names="kinetic_id", level="well", inplace=True)
    return df_kinetics


def concatenate_kinetics(*, run_list, runinfo_df, fp_experiment, cut_cycle=False):
    """Concatenate the DataFrames for different rounds.

    Parameters
    ----------
    run_list : list
        List of all run_ids to be analysed 
    runinfo_df : pandas.Dataframe
        Index is the number of rounds, necessary columns are run_id-VK, result_file-VK, run_id-HK, result_file-HK, run_id-Assay, result_file(s)-Assay.
    fp_experiment : pathlib.Path
        The path to the top-level folder containing all run folders.
    cut_cycle : optional, int
        All measurements must have the same number of cycles. If not, the cut_cycle determines where the measurements are cut.
    concentration_factor : int
        Default 1, can be adapted if different between runs
    """ 
    dfs_kinetics = [
        read_kinetics_run(runinfo_df=runinfo_df, fp_experiment=fp_experiment, round_no=round_no)
        for round_no, _ in enumerate(run_list)
    ]
    df_kinetics = pandas.concat(dfs_kinetics)
    if cut_cycle:
        df_kinetics = df_kinetics[df_kinetics.index.get_level_values("cycle").isin(numpy.arange(cut_cycle+1))]
    return df_kinetics
