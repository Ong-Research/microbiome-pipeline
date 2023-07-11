
configfile: "config/config.yaml"
targets = [config["metadata_file"], config["negcontrol_file"]]


rule all:
  input: targets

rule prepare_metadata:
  """Prepares a metadata tibble based on previous data"""
  input:
    config["samples_file"]
  output:
    config["metadata_file"]
  script:
    """./scripts/init/prepare_metadata.R"""


rule prepare_negcontrols:
  """Prepares a tibble with the negative controls that correspond to each sample"""
  input:
    samples = config["samples_file"],
    meta = config["metadata_file"]
  params:
    batch_2018 = "data/negcontrols_meta/2018/Nextseq_04202018_mapping_WISC.txt",
    batch_2020 = "data/negcontrols_meta/2020/Nextseq_190826_mapping_WISC.txt",
    batch_2022 = "data/negcontrols_meta/2022/WISC_220214_seqmap.txt"
  output:
    config["negcontrol_file"]
  script:
    """./scripts/init/prepare_negcontrols.R"""
