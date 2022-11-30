#!/usr/bin/env python
# coding=utf-8

# This module attempts to gather and represent info from _coverage_metrics.csv files. It relies on the official standard: 
# https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/CoverageMetricsReport_fDG.htm
# Still in development. 

import re
import logging

from multiqc.utils.config import read_count_multiplier, base_count_multiplier, read_count_prefix, base_count_prefix

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table
from multiqc.utils.util_functions import write_data_file


# Initialise the logger
log = logging.getLogger(__name__)


NAMESPACE = "DRAGEN coverage"


class DragenCoverageMetrics(BaseMultiqcModule):
    def add_coverage_metrics(self):

        # Holders for the table's data parameter.
        sample_by_phenotype = {}
        phenotype_by_metric = {}

        for f in self.find_log_files("dragen/coverage_metrics"):
        
            sample, phenotype, data = parse_coverage_metrics(f)
            
            if not sample: continue

            sample = self.clean_s_name(sample, f) 
            clean_sp(sample, phenotype)
            
            self.add_data_source(f, section = "stats") 
            
            if sample not in phenotype_by_metric: phenotype_by_metric[sample] = {}
            
            combo_name = sample + "_" + phenotype
            if combo_name not in sample_by_phenotype: sample_by_phenotype[combo_name] = {}
            
            for metric in data:
                value = data[metric]
                sample_by_phenotype[combo_name][metric] = value
                
                combo_name2 = phenotype + " " + metric
                if combo_name2 not in phenotype_by_metric[sample]: phenotype_by_metric[sample][combo_name2] = {}
                phenotype_by_metric[sample][combo_name2] = value
                

        # Filter to strip out ignored sample names.
        sample_by_phenotype = self.ignore_samples(sample_by_phenotype)

        # Write data to file
        self.write_data_file(sample_by_phenotype, "dragen_cov_metrics")

        self.general_stats_addcols(sample_by_phenotype, table_headers_sample_by_phenotype, namespace = "DRAGEN coverage")


        """
        The most interesting parameter of the following method is "plot", which is basically just a html code block.
        A perfect place for js injection. Could be useful to transfer data: "<script>DRAGEN_data={some_var}</script>".format(some_var = json.dumps(some_var)) 
        CSS/JS could be also imported, see fastqc.py for details.      
        """

        self.add_section(
            name="Coverage metrics",
            anchor="dragen-cov-metrics",
            description="""
            Coverage metrics over a region (where the region can be a target region,
            a QC coverage region, or the whole genome). Press the `Help` button for details.
            """,
            helptext="""
            The following criteria are used when calculating coverage:

            * Duplicate reads and clipped bases are ignored.
            * Only reads with `MAPQ` > `min MAPQ` and bases with `BQ` > `min BQ` are considered

            Considering only bases usable for variant calling, _i.e._ excluding:

            1. Clipped bases
            2. Bases in duplicate reads
            3. Reads with `MAPQ` < `min MAPQ` (default `20`)
            4. Bases with `BQ` < `min BQ` (default `10`)
            5. Reads with `MAPQ` = `0` (multimappers)
            6. Overlapping mates are double-counted
            """,
            plot = table.plot(sample_by_phenotype, table_headers_sample_by_phenotype, pconfig = {"namespace": NAMESPACE}) +
                   table.plot(phenotype_by_metric, table_headers_phenotype_by_metric, pconfig = {"namespace": NAMESPACE}),

        )
        return sample_by_phenotype.keys()
    
    
    
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The following block aims at automating the parsing of coverage metrics. 
It also gives some direct control over the final html table view.

"cov_pats" dictionary:

    Keys:   regexes to match a single metric or a group of metrics.

    Values: "gen_table" returns a string corresponding to the official metric's name.
            "own_table"
            "configs" is a dict with configs for html table. Most values are based on previous version of this file and /dragen/utils.py  
    
    
Available configs: https://github.com/ewels/MultiQC/blob/master/docs/plots.md#creating-a-table


Notes: 

- All values of "metric" and "configs" must be lambda x: or normal functions for consistency and parsing convinience. 
  They must return some valid (str for "metric" and link above for "configs") parameter: lambda **x: par 1
  If metrics have two values (only "Aligned bases/reads in region" metrics for now), then they must return: lambda **x: [par 1, par 2].
  This is needed because we need to split these values in two columns and modify them separately.

- Every pattern must contain exactly one group, because it is a parameter for the lambda functions.
  () at the end of some patterns (e.g. matches only single metric) is used to match an empty string just for group to exist.

- All functions have only one argument **x. How it is handled depends on the lambdas' expressions.

- The "cov_pats" version below is the most stable at the moment but old and will look  somewhat different in future.

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

cov_pats = {

"Aligned (bases|reads)$": {
    "metric":          lambda x: "Aligned " + x + ".",
    "configs":{
        "description": lambda x: "Total number of aligned " + x,
        "title":       lambda x: read_count_prefix + " Aln " + x if x == "reads" else base_count_prefix + " Aln " + x,
        "scale":       lambda x: "RdYlGn" if x == "bases" else "Greys",
        "format":      lambda x: "{:,.0f}",
        "modify":      lambda x: (lambda y: y * read_count_multiplier) if x == "reads" else (lambda y: y * base_count_multiplier),
        "hidden":      lambda x: True,
    }
},
"Aligned (bases|reads) in .+": {
    "metric":           lambda x: [lambda x: "Aligned " + x + " in region.", lambda x: "Percentage of aligned " + x + " in region."],
    "configs":{
        "description": [lambda x: "Number of uniquely mapped bases to region." if x == "bases" else "Number of uniquely mapped reads to region. When region is the target BED, this metric is equivalent to and replaces Capture Specificity based on target region.",
                        lambda x: "Percentage relative to the number of uniquely mapped bases to the genome." if x == "bases" else "Percentage relative to the number of uniquely mapped reads to the genome. When region is the target BED, this metric is equivalent to and replaces Capture Specificity based on target region."],
        "title":       [lambda x: "{} Aln ".format(read_count_prefix if x == "reads" else base_count_prefix) + x + " on trg",
                        lambda x: "{} Aln ".format(read_count_prefix if x == "reads" else base_count_prefix) + x + " on trg"],
        "scale":       [lambda x: "RdGy" if x == "bases" else "RdYlGn", lambda x: "Greys" if x == "bases" else "Greens"],
        "format":      [lambda x: "{:,.0f}", lambda x: "{:,.0f}"],
        "suffix":      [lambda x: "", lambda x: " %"],
        "max":         [lambda x: 100, lambda x: 100],
        "modify":      [lambda x: (lambda y: y * read_count_multiplier) if x == "reads" else (lambda y: y * base_count_multiplier),
                        lambda x: (lambda y: y * read_count_multiplier) if x == "reads" else (lambda y: y * base_count_multiplier)],
    }
},
"Average alignment coverage over .+()": {
    "metric":          lambda x: "Average alignment coverage over region.",
    "configs":{
        "description": lambda x: "Number of uniquely mapped bases to region divided by the number of sites in region.",
        "title":       lambda x: "Depth",
        "scale":       lambda x: "Blues",
        "suffix":      lambda x: " x",
        "format":      lambda x: "{:,.1f}",
    }
},
"Uniformity of coverage \((.+)\) over .+": {
    "metric":          lambda x: "Uniformity of coverage " + x + " over region",
    "configs":{
        "description": lambda x: "Percentage of sites with coverage greater than % of the mean coverage in region.",
        "title":       lambda x: re.search("PCT(.+)", x, re.IGNORECASE).group(1).replace(" ", ""),
        "scale":       lambda x: "RdBu",
        "suffix":      lambda x: " %",
        "format":      lambda x: "{:,.1f}",
    }
},
"PCT of .+ with coverage (.+)": {
    "metric":          lambda x: "PCT of region with coverage region " + x,
    "configs":{
        "description": lambda x: "Percentage of sites in region with at least " + re.search("(\d+)x:", x).group(1) + " but less than " + re.search(":\s*(.+)\)", x).group(1) + " coverage.",
        "title":       lambda x: x.replace(" ", ""),
        "scale":       lambda x: "RdBu",
        "suffix":      lambda x: " %",
        "format":      lambda x: "{:,.1f}",
    }
},
"Average (chr X|chr Y|mitochondrial|autosomal) coverage over .+": {
    "metric":          lambda x: "Average " + x + " coverage over region.",
    "configs":{
        "description": lambda x: "Total number of bases that aligned to the intersection of chromosome Y with region divided by the total number of loci in the intersection of chromosome Y with region. If there is no chromosome Y in the reference genome or the region does not intersect chromosome Y, this metric shows as NA.",
        "title":       lambda x: "Avg " + re.search("(X|Y|mit|aut)", x, re.IGNORECASE).group(1) + " cov",
        "scale":       lambda x: "Blues",
        "suffix":      lambda x: " x",
        "format":      lambda x: "{:,.1f}",
    }
},
"Median autosomal coverage over .+()": {
    "metric":          lambda x: "Median autosomal coverage over region.",
    "configs":{
        "description": lambda x: "Median alignment coverage over the autosomal loci in region. If there is no autosome in the reference genome or the region does not intersect autosomes, this metric shows as NA.",
        "title":       lambda x: "Med aut cov",
        "scale":       lambda x: "Blues",
        "suffix":      lambda x: " x",
        "format":      lambda x: "{:,.1f}",
    }
},
"Mean/Median autosomal coverage ratio over .+()": {
    "metric":          lambda x: "Mean/Median autosomal coverage ratio over region.",
    "configs":{
        "description": lambda x: "Mean autosomal coverage in region divided by the median autosomal coverage in region. If there is no autosome in the reference genome or the region does not intersect autosomes, this metric shows as NA.",
        "title":       lambda x: "Mean/Med aut cov",
    }
},
"XAvgCov/YAvgCov ratio over .+()": {
    "metric":          lambda x: "XAvgCov/YAvgCov ratio over region.",
    "configs":{
        "description": lambda x: "Average chromosome X alignment coverage in region divided by the average chromosome Y alignment coverage in region. If there is no chromosome X or chromosome Y in the reference genome or the region does not intersect chromosome X or Y, this metric shows as NA.",
        "title":       lambda x: "XAvgCov/YAvgCov",
    }
},
"(X|Y)AvgCov/AutosomalAvgCov ratio over .+": {
    "metric":          lambda x: x + "AvgCov/AutosomalAvgCov ratio over region.",
    "configs":{
        "description": lambda x: "Not given in the official standard. Please fill this in if possible",
        "title":       lambda x: x + "AvgCov/AutAvgCov",
    }
}}

# Add common configs in a loop.
for metric_pattern in cov_pats:
    temp = cov_pats[metric_pattern]["configs"]
    if isinstance(temp["description"], list):
        temp["min"] = [lambda x: 0, lambda x: 0]
        temp["hidden"] = [lambda x: True, lambda x: True]
    else:
        temp["min"] = lambda x: 0
        temp["hidden"] = lambda x: False




table_headers_sample_by_phenotype = {}
table_headers_phenotype_by_metric = {}

def parse_coverage_metrics(f):

    """
    Some possible file names:
    
    T_SRR7890936_50pc.wgs_coverage_metrics_normal.csv
    T_SRR7890936_50pc.wgs_coverage_metrics_tumor.csv
    T_SRR7890936_50pc.target_bed_coverage_metrics.csv
    SAMPLE_MONO_0_dragen.qc-coverage-region-1_coverage_metrics.csv
    SAMPLE_MONO_1_dragen.wgs_coverage_metrics.csv
    """

    r"""
    The previous search pattern from \multiqc\utils\search_patterns.yaml: ".*\.(wgs|target_bed)_coverage_metrics_?(tumor|normal)?\.csv" was not valid.
    The new search regex has the most general structure : ".*_coverage_metrics.*\.csv"
    File name is checked below, because there is no guarantee that it would be as in the official standard: _coverage_metrics.csv
    All unusual names will be reported, after that the following block can be modified. 
    """

    m = re.search(r"(.+)\.(.+)_coverage_metrics(.*).csv", f["fn"])
    if not m:
    
        # Check for the absent point.
        m = re.search(r"(.+)()_coverage_metrics(.*).csv", f["fn"])
        if not m:
            tempStr = " " * 88
            log.debug("Only coverage metric files of the following form are supported:\n" + tempStr +
                      "+_coverage_metrics*.csv\n" + tempStr +
                      "where + is some non-empty string and * may be empty.\n" + tempStr +
                      "The following file is not valid:\n" + tempStr + 
                      "%s root: %s", f['fn'], f['root'])
            return 0, 0, 0, 0

    sample, phenotype = m.group(1), m.group(2)
    if m.group(3): phenotype += '_' + m.group(3) # Both parts are just concatenated together for simplicity. Some other handling could be more appropriate.

    sample, phenotype = clean_sp(sample, phenotype)
    #--------------------------------------------------------------#
    global table_headers_sample_by_phenotype
    global table_headers_phenotype_by_metric
    data = {}

    for line in f["f"].splitlines():
    
        # Check the general structure of line.
        m = re.search("COVERAGE SUMMARY,,([^,]+),([^,]+),?(.*)", line, re.IGNORECASE)

        # If line is fine then extract the necessary fields.
        if m:
            metric, value1, value2 = m.group(1), m.group(2), m.group(3)
            
        # Otherwise report and go to the next line.
        else:
            log.debug("An extraordinary line was found in %s: %s", f['fn'], line)
            continue
        
        
        # Check values. If not int or float report to the log. 
        # If some unusual values are found later, then this block can be extended to handle them.
        try:
            value1 = int(value1)
        except ValueError:
            try:
                value1 = float(value1)
            except ValueError:
                log.debug("An extraordinary value (non int/float) was found in %s: %s, %s, %s", f['fn'], metric, value1, value2)
                
        if value2:
            try:
                value2 = float(value2)
            except ValueError:
                try:
                    value2 = int(value2)
                except ValueError:
                    log.debug("An extraordinary value (non int/float) was found in %s: %s, %s, %s", f['fn'], metric, value1, value2)
                    

        """
        "metric" is checked in the following.
        """

        done = 0

        for metric_pattern in cov_pats:

            metric_match = re.search(metric_pattern, metric, re.IGNORECASE)           
            if metric_match:

                metric_match = metric_match.group(1)
                metric = cov_pats[metric_pattern]["metric"](metric_match)
                metric_pct = ''

                if isinstance(metric, list):
                    metric_pct = metric[1](metric_match)
                    metric = metric[0](metric_match)

                if metric not in table_headers_sample_by_phenotype:
                
                    table_headers_sample_by_phenotype[metric] = {}
                    if metric_pct: table_headers_sample_by_phenotype[metric_pct] = {}

                    for header, head_value in cov_pats[metric_pattern]["configs"].items():
                    
                        if isinstance(head_value, list):
                            table_headers_sample_by_phenotype[metric][header] = head_value[0](metric_match)
                            if metric_pct:  table_headers_sample_by_phenotype[metric_pct][header] = head_value[1](metric_match)
                        else:
                            table_headers_sample_by_phenotype[metric][header] = head_value(metric_match)
                  
                metric_pct = ''
                if isinstance(metric, list): metric_pct = metric[1](metric_match)

                metric2 = phenotype + " " + metric
                if metric2 not in table_headers_phenotype_by_metric:
                
                    table_headers_phenotype_by_metric[metric2] = {}
                    if metric_pct: table_headers_phenotype_by_metric[metric_pct] = {}

                    for header, head_value in cov_pats[metric_pattern]["configs"].items():
                    
                        if isinstance(head_value, list):
                            table_headers_phenotype_by_metric[metric2][header] = head_value[0](metric_match)
                            if metric_pct:  table_headers_phenotype_by_metric[metric_pct][header] = head_value[1](metric_match)
                            table_headers_phenotype_by_metric[metric2]["hidden"] = False
                            #if metric_pct:
                            if metric_pct:  table_headers_phenotype_by_metric[metric_pct][header] = head_value[1](metric_match)
                            table_headers_phenotype_by_metric[metric2]["title"] =table_headers_sample_by_phenotype[metric]['title']+ phenotype

                        else:
                            table_headers_phenotype_by_metric[metric2][header] = head_value(metric_match)
                            table_headers_phenotype_by_metric[metric2]["title"] =table_headers_sample_by_phenotype[metric]['title']+ ' '+phenotype
                            table_headers_phenotype_by_metric[metric2]["hidden"] = False
                
                data[metric] = value1
                if value2: data[metric_pct] = value2
                
                done = 1
                break
                
        # Final case: no pattern has matched the "metric".
        if not done:
            log.debug("An extraordinary metric was found in %s: %s", f['fn'], metric)
            continue
            
    return sample, phenotype, data

# Some silly cleaning for now. 
# Probably better to write some common dragen cleaner for consistency across dragen modules.
# One could also extend the base_module.clean_s_name. Not sure though if touching the global cleaner would be a good idea.
def clean_sp(sample, phenotype):
    
    sample = sample.replace("-", " ")
    sample = sample.replace("_", " ")
    phenotype = phenotype.replace("-", " ")
    phenotype = phenotype.replace("_", " ")
    
    phenotype = phenotype.replace("region", "")
    phenotype = phenotype.replace("coverage", "")
    sample = sample.replace("dragen", "")
    sample = sample.replace("SAMPLE", "")
    
    return sample, phenotype
