'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This module gathers CNV metrics data and prepares it for the output report.

It relies on the following official source:
https://support-docs.illumina.com/SW/DRAGEN_v310/Content/SW/DRAGEN/CNVMetrics.htm
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""'''

import re
import logging
from collections import defaultdict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table
from multiqc.utils.util_functions import write_data_file

"""
The common dragen cleaner can be written later, for now it is placed locally,
but is disabled to provide the compatibility with other dragen modules.
from .utils import clean_sample
"""


NAMESPACE = "DRAGEN CNV"


"""
The SINGLE_HEADER is used to define the common settings for all headers:
https://github.com/ewels/MultiQC/blob/master/docs/plots.md#creating-a-table
"""
SINGLE_HEADER = {
    "namespace": NAMESPACE,        # Name for grouping. Prepends desc and is in Config Columns modal
    "title": "",                   # Short title, table column title
    "description": "",             # Longer description, goes in mouse hover text
    "max": None,                   # Minimum value in range, for bar / colour coding
    "min": 0,                      # Maximum value in range, for bar / colour coding
    "ceiling": None,               # Maximum value for automatic bar limit
    "floor": None,                 # Minimum value for automatic bar limit
    "minRange": None,              # Minimum range for automatic bar
    "scale": "GnBu",               # Colour scale for colour coding. False to disable.
    "bgcols": {},                  # Dict with values: background colours for categorical data.
    "colour": "255, 150, 30",      # Colour for column grouping
    "suffix": "",                  # Suffix for value (eg. "%")
    "format": "{:,.1f}",           # Value format string - default 1 decimal place
    "cond_formatting_rules": {},   # Rules for conditional formatting table cell values.
    "cond_formatting_colours": [], # Styles for conditional formatting of table cell values
    "shared_key": None,            # See the link for description
    "modify": None,                # Lambda function to modify values, special case, see below
    "hidden": True,                # Set to True to hide the column on page load
    "hidden2": False,              # For non-general plots in the CNV section.
}
# The EXTRA_HEADER contains the logical keys of the SINGLE_HEADER.
# These are just extensions, but can be used as if they were a part of the SINGLE_HEADER.
# These are not set by default. Only the SINGLE_HEADER is used as a basis for a header.
EXTRA_HEADER = {
    "hidden2": False,              # For non-general plots in the CNV section.
}


'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
You can set any possible real/virtual configuration defined in the SINGLE_HEADER.

Some rules:

- "modify" must be (lambda x: x if isinstance(x, str) else expression)
  Or generally speaking your function must check for string and return it.
  If not string (eg int or float), then do something with x and return the result.

- if you wish to add new "CNV SUMMARY metrics" please take a look at the "make_metric_id".
'''

V2 = "(second value)"  # Used to create a unique ID if a metric has two values.

SEX_GENOTYPER_SECTION = "SEX GENOTYPER" # "Sex Genotyper Metrics" section.
CNV_SUMMARY_SECTION   = "CNV SUMMARY"   # "CNV Summary   Metrics" section.

# Used to make a number "narrower" to not overlap with adjacent right columns.
_modify = lambda x: x if isinstance(x, str) else x * 0.000001

METRICS = {
    SEX_GENOTYPER_SECTION:{
        "karyotype":{
            "title": "Kar", "min": None, "colour": "0, 255, 0", "hidden": False,
            "description": "Estimated sex karyotype for the sample.",
            "cond_formatting_rules":{
                "red": [{"s_ne": "XX"},{"s_ne": "YX"},{"s_ne": "XY"}],
                "green": [{"s_eq": "XX"},{"s_eq": "YX"},{"s_eq": "XY"}],
            },
            "cond_formatting_colours":[
                {"red": "#FF0000"},
                {"green": "#00FF00"},
            ],
        },
        "confidence":{
            "title": "KarConf", "colour": "0, 255, 0", "hidden": False,
            "suffix": " %", "format": "{:,.3f}", "max": 1.0,
            "cond_formatting_rules":{
                "red": [{"lt": 0},{"gt": 1}],
            },
            "cond_formatting_colours":[
                {"red": "#FF0000"},
            ],
            "description": "Confidence metric ranging from 0.0 to 1.0. " +
                           "If the sample sex is specified, this metric is 0.0.",
        },
    },
    CNV_SUMMARY_SECTION:{
        "bases in reference genome":{
            "title": "Mb bases", "scale": "RdYlGn", "colour": "0, 0, 255", "format": "{:,.0f}",
            "description": "Bases in reference genome in use.", "modify": _modify,
        },
        "average alignment coverage over genome":{
            "title": "AvgCov", "suffix": " %",
            "description": "Average alignment coverage over genome.",
        },
        "number of alignment records":{
            "title": "M Aln records", "format": "{:,.0f}", "modify": _modify,
            "description": "Number of alignment records processed.",
        },
        "coverage uniformity":{
            "title": "CovUnif", "suffix": " %", "format": "{:,.2f}",
            "description": "",
        },
        "number of target intervals":{
            "title": "Trg intervals", "format": "{:,.0f}",
            "description": "Number of target intervals.",
        },
        "number of segments":{
            "title": "Segments", "format": "{:,.0f}",
            "description": "Number of segments.",
        },


        "number of filtered records (total)":{
            "title": "M Total", "format": "{:,.2f}", "modify": _modify,
            "description": "Number of filtered records (total).",
        },
        "number of filtered records (total)" + V2:{
            "title": "Total%", "suffix": " %",
            "description": "Percentage of filtered records (total).",
        },
        "number of filtered records (duplicates)":{
            "title": "M Duplicates", "format": "{:,.2f}",  "modify": _modify,
            "description": "Number of filtered records (due to duplicates).",
        },
        "number of filtered records (duplicates)" + V2:{
            "title": "Duplicates%", "suffix": " %",
            "description": "Percentage of filtered records (due to duplicates).",
        },
        "number of filtered records (mapq)":{
            "title": "M MAPQ", "format": "{:,.2f}", "modify": _modify,
            "description": "Number of filtered records (due to MAPQ).",
        },
        "number of filtered records (mapq)" + V2:{
            "title": "MAPQ%", "suffix": " %",
            "description": "Percentage of filtered records (due to MAPQ).",
        },
        "number of filtered records (unmapped)":{
            "title": "M Unmapped", "format": "{:,.2f}", "modify": _modify,
            "description": "Number of filtered records (due to being unmapped).",
        },
        "number of filtered records (unmapped)" + V2:{
            "title": "Unmapped%", "suffix": " %",
            "description": "Percentage of filtered records (due to being unmapped).",
        },


        "number of amplifications":{
            "title": "Amplifs", "format": "{:,.0f}",
            "description": "Number of amplifications.",
        },
        "number of passing amplifications":{
            "title": "Amplifs PASS", "format": "{:,.0f}",
            "description": "Number of PASS amplifications.",
        },
        "number of passing amplifications" + V2:{
            "title": "Amplifs PASS%", "suffix": " %",
            "description": "Percentage of PASS amplifications.",
        },


        "number of deletions":{
            "title": "Deletions", "format": "{:,.0f}",
            "description": "Number of deletions.",
        },
        "number of passing deletions":{
            "title": "Deletions PASS", "format": "{:,.0f}",
            "description": "Number of PASS deletions.",
        },
        "number of passing deletions" + V2:{
            "title": "Deletions PASS%", "suffix": " %",
            "description": "Percentage of PASS deletions.",
        },


        "number of de novo calls":{
            "title": "de Novo", "format": "{:,.0f}",
        },
        "number of passing de novo calls":{
            "title": "de Novo PASS", "format": "{:,.0f}",
            "description": "Number of passing de novo calls.",
        },
        "number of passing de novo calls" + V2:{
            "title": "de Novo PASS%", "suffix": " %",
            "description": "Percentage of passing de novo calls.",
        },
    },
}



"""
The dragen.py provides a common interface for all dragen modules.
To avoid cluttering up the global dragen space, only one method is defined.
Other methods can be added as well to provide modularity.
Collected CNV data can be also transferred to other dragen modules as an attribute.
"""
cnv_data = {}
cnv_headers = {}

class DragenCNVMetrics(BaseMultiqcModule):
    def add_cnv_metrics(self):
        """Called from dragen.py. Does different abstract tasks."""

        global cnv_data

        for file in self.find_log_files("dragen/cnv_metrics"):
            if parse(file): self.add_data_source(file, section = "stats")

        if not cnv_data: return set()

        cnv_data = self.ignore_samples(cnv_data)
        '''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        In theory some combinations of sample names, configurations and implementation aspects
        of both cleaners may produce duplicates. Cleaned names are checked for duplicates,
        if presented then their original names are used instead. The issue will be reported.
        """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""'''
        check_cleaned_names(
        {s_name: clean_sample(self.clean_s_name(s_name)) for s_name in cnv_data.keys()})

        # Write data into the general table.
        self.general_stats_addcols(cnv_data, make_general_headers(), namespace = NAMESPACE)

        # Write data to file.
        self.write_data_file(cnv_data, "dragen_cnv_metrics")

        # Report all found errors and warnings.
        make_log_report()

        self.add_section(
            name = "CNV metrics",
            anchor = "dragen-cnv-metrics",
            description = """
            DRAGEN CNV outputs metrics in CSV format.
            The output follows the general convention for QC metrics reporting in DRAGEN.
            The CNV metrics are output to a file with a *.cnv_metrics.csv file extension.
            The following list summarizes the metrics that are output from a CNV run.
            Press the `Help` button for details.
            """,
            helptext = """
            Sex Genotyper Metrics:

            * Estimated sex karyotype for the sample, along with a confidence metric ranging from 0.0 to 1.0.
            If the sample sex is specified, this metric is 0.0.
            * If using a panel of normals, all panel samples are also reported.

            CNV Summary Metrics:

            * Bases in reference genome in use.
            * Average alignment coverage over genome.
            * Number of alignment records processed.
            * Number of filtered records (total).
            * Number of filtered records (due to duplicates).
            * Number of filtered records (due to MAPQ).
            * Number of filtered records (due to being unmapped).
            * Number of target intervals.
            * Number of normal samples.
            * Number of segments.
            * Number of amplifications.
            * Number of deletions.
            * Number of PASS amplifications.
            * Number of PASS deletions.
            """,
            plot = make_own_plot(),
        )
        return cnv_data.keys()



"""
The log_data is the common holder for errors and warnings,
which are stored in groups to make the output log file more readable.
errors:
- duplicate_file_names, which theoretically can be found in subdirectories.
- duplicate_clean_sample_names is a special case for catching bugs associated with cleaning.
  Details are in the DragenCNVMetrics.
- invalid_file_lines, which do not conform to:
    SEX GENOTYPER,,<sample name>,<sex karyotype>,<confidence> or
    CNV SUMMARY,,<metric>,<value1>  or
    CNV SUMMARY,,<metric>,<value1>,<value2>

warnings:
- unknown_metrics are not presented in the METRICS.
- unusual_values are those except for int and float in CNV SUMMARY section.
                -//-              XY-system/float in SEX GENOTYPER section.
- unknown_sections can be potentially implemented.
  For now "invalid_file_lines" is used instead.
"""
log = logging.getLogger(__name__)

log_data = {
    "duplicate_file_names": defaultdict(list),
    "duplicate_clean_sample_names": {},
    "invalid_file_lines": defaultdict(lambda: defaultdict(list)),
    "unknown_metrics": defaultdict(lambda: defaultdict(list)),
    "unusual_values": defaultdict(lambda: defaultdict(dict)),
    "unknown_sections": defaultdict(dict),
}


def parse(file_handler):
    """Parser for cnv_metrics.csv files.
    Input:  file_handler.
    Output: 0 as failure or 1 as success in extracting of metrics."""

    root_name = file_handler["root"]
    file_name = file_handler["fn"]

    # Official structure of files: *.cnv_metrics.csv
    # Exepted structure: <arbitrary prefix>.cnv_metrics<arbitrary suffix>.csv
    file_match = re.search(r"(.*)\.cnv_metrics(.*).csv", file_name)

    sample, arb_suffix = file_match.group(1), file_match.group(2)
    if not sample: sample = "Unknown sample"
    if arb_suffix: sample += arb_suffix

    if sample in cnv_data:
        log_data["duplicate_file_names"][root_name].append(file_name)
        return 0

    success = 0
    data = {}
    extra_count = -1 # Special counter to handle the case when multiple rows are presented.

    for line in file_handler["f"].splitlines():
        # Check the general structure.
        line_match = re.search("([^,]+),,([^,]+),([^,]+),?([^,]*)", line)

        # If line is fine then extract the necessary fields.
        if line_match:
            section, metric, value1, value2 = line_match.group(1), line_match.group(2),\
                                              line_match.group(3), line_match.group(4)

            if re.search("^\s*CNV\s*SUMMARY\s*$", section, re.IGNORECASE):

                metric_id = make_metric_id(metric)

                """ Check the extracted values.
                Type conversion is only performed to int and float.
                Other types of data are left unchanged as strings."""
                try:
                    value1 = int(value1)
                except ValueError:
                    try:
                        value1 = float(value1)
                    except ValueError:
                        log_data["unusual_values"][root_name][file_name][metric] = value1
                data[metric_id] = value1

                if metric_id not in cnv_headers:
                    cnv_headers[metric_id] = make_header(metric_id, CNV_SUMMARY_SECTION)

                if value2:
                    try:
                        value2 = float(value2)
                    except ValueError:
                        try:
                            value2 = int(value2)
                        except ValueError:
                            log_data["unusual_values"][root_name][file_name][metric + V2] = value2
                    metric_id += V2
                    data[metric_id] = value2
                    if metric_id not in cnv_headers:
                        cnv_headers[metric_id] = make_header(metric_id, CNV_SUMMARY_SECTION)

            elif re.search("^\s*SEX\s*GENOTYPER\s*$", section, re.IGNORECASE):
                """
                SEX GENOTYPER,,SAMPLE_MONO_1,XX,0.9959
                The code supports only the metric of this type.

                The metric(3rd col) itself is just a subset of a sample name,
                hence can not be used as a special common metric name.
                """
                # In case the metric(sample name) is also produced by some other module.
                metric_id = metric + "_CNV_metrics"
                extra_count += 1
                if extra_count: metric_id += str(extra_count)

                special_case = 0
                # Does it at least contain X or Y?
                if re.search("[XY]", value1, re.IGNORECASE):
                    special_case = 1
                    metric_id = "Karyotype"
                    cnv_headers[metric_id] = make_header("karyotype", SEX_GENOTYPER_SECTION)
                else:
                    cnv_headers[metric_id] = make_header(metric_id, SEX_GENOTYPER_SECTION)

                data[metric_id] = value1

                # The regex handles the general line structure,
                # which does not necessarily has 2 values.
                if not value2: continue

                try:
                    value2 = float(value2)
                except ValueError:
                    try:
                        value2 = int(value2)
                    except ValueError:
                        pass

                if special_case: metric_id = "Confidence"
                else:
                    metric_id += V2
                    if extra_count: metric_id += str(extra_count)

                if metric_id not in cnv_headers:
                    if special_case:
                        cnv_headers[metric_id] = make_header("confidence", SEX_GENOTYPER_SECTION)
                    else:
                        cnv_headers[metric_id] = make_header(metric_id, CNV_SUMMARY_SECTION)
                data[metric_id] = value2
            else:
                #log_data["unknown_sections"][root_name][file_name] = line
                log_data["invalid_file_lines"][root_name][file_name].append(line)
                continue

            success = 1

        # Otherwise check if line is empty. If not then report it and go to the next line.
        else:
            if not re.search("^\s*$", line):
                log_data["invalid_file_lines"][root_name][file_name].append(line)
            continue

    if success: cnv_data[sample] = data
    return success


def make_metric_id(metric):
    """Improves consistency to a certain degree."""
    # Metric-specific peculiarities can be handled here.

    # Single backslashes are escaped if presented.
    metric_id = re.sub("\s+", " ", metric).strip()
    metric_id = metric_id.lower()

    return metric_id


def make_header(metric, section):
    """Creates a single header for a given metric."""
    configs = {}

    for key, val in SINGLE_HEADER.items():
        if val is not None: configs[key] = val

    if section in METRICS:
        if metric in METRICS[section]:
            configs.update(METRICS[section][metric])

    if not("title" in configs and configs["title"]):
        configs["title"] = metric

    return configs


def make_general_headers():
    """Prepares headers for the general section."""
    headers = {}

    for metric in cnv_headers:
        headers[metric] = {
            key: val for key, val in cnv_headers[metric].items() if key not in EXTRA_HEADER
        }
    return headers

def make_own_plot():
    """Creates a plot for the CNV section."""

    for header in EXTRA_HEADER:
        for metric in cnv_headers:
            if header in cnv_headers[metric]:
                if header == "hidden2":
                    cnv_headers[metric]["hidden"] = cnv_headers[metric]["hidden2"]
                del cnv_headers[metric][header]

    return table.plot(cnv_data, cnv_headers, {"namespace": NAMESPACE})


##-----------------------------------------------------------------------------##
##                              Various functions                              ##
##-----------------------------------------------------------------------------##

def check_cleaned_names(sample_names):
    """Check cleaned sample names for duplicates."""

    temp = defaultdict(list)
    for orig_s_name in sample_names: temp[sample_names[orig_s_name]].append(orig_s_name)

    for clean_s_name in temp:
        if len(temp[clean_s_name]) > 1:
            log_data["duplicate_clean_sample_names"][clean_s_name] = temp[clean_s_name]
        else:
            orig_s_name = temp[clean_s_name][0]
            if clean_s_name != orig_s_name: cnv_data[clean_s_name] = cnv_data.pop(orig_s_name)


def make_log_report():
    """Creates a readable and informative log output."""

    if log_data["duplicate_file_names"]:
        log_message = "\n\nThe following duplicates were found:\n"

        for root in log_data["duplicate_file_names"]:
            log_message += "  " + root + ":\n"
            for file in log_data["duplicate_file_names"][root]:
                log_message += "    " + file + "\n"

        log.error(log_message + "\n")

    if log_data["duplicate_clean_sample_names"]:
        log_message = ("\n\nThe following sample names are left as they are " +
                       "because their cleaned versions are duplicates.\n")

        for clean_s_name in log_data["duplicate_clean_sample_names"]:
            log_message += "  " + clean_s_name + " is formed from:\n"
            for orig_s_name in log_data["duplicate_clean_sample_names"][clean_s_name]:
                log_message += "    " + orig_s_name + "\n"

        log.error(log_message + "\n")

    if log_data["invalid_file_lines"]:
        log_message = ("\n\nThe lines in files must be:\n" +
                       "CNV SUMMARY,,<metric>,<value1> or\n" +
                       "CNV SUMMARY,,<metric>,<value1>,<value2> or\n" +
                       "SEX GENOTYPER,,<metric>,<value1>,<value2>\n\n"  +
                       "The following files contain invalid lines:\n" )

        for root in log_data["invalid_file_lines"]:
            log_message += "  " + root + ":\n"
            for file in log_data["invalid_file_lines"][root]:
                log_message += "    " + file + ":\n"
                for line in log_data["invalid_file_lines"][root][file]:
                    log_message += "      " + line + "\n"

        log.error(log_message + "\n")

    if log_data["unknown_metrics"]:
        log_message = "\n\nThe following files contain unknown metrics:\n"

        for root in log_data["unknown_metrics"]:
            log_message += "  " + root + ":\n"
            for file in log_data["unknown_metrics"][root]:
                log_message += "    " + file + ":\n"
                for metric in log_data["unknown_metrics"][root][file]:
                    log_message += "      " + metric + "\n"

        log.warning(log_message + "\n")

    if log_data["unusual_values"]:
        log_message = ("\n\nAll metrics' values except int, float and NA are non-standard.\n"  +
                       "The following files contain non-standard values:\n" )

        for root in log_data["unusual_values"]:
            log_message += "  " + root + ":\n"
            for file in log_data["unusual_values"][root]:
                log_message += "    " + file + ":\n"
                for metric in log_data["unusual_values"][root][file]:
                    log_message += ("      " + metric + " = " +
                                    log_data["unusual_values"][root][file][metric] + "\n")

        log.warning(log_message + "\n")

# This cleaner will be exported into the ./dragen/utils.py later.
# Turned off for compatibility with other dragen modules.
def clean_sample(sample):
    """Cleans sample names to make the output tables more attractive."""
    return sample

    sample = re.sub("dragen|sample", "", sample, flags = re.I)
    sample = re.sub(r"-|_|\.", " ", sample)
    sample = sample.strip()

    return sample
