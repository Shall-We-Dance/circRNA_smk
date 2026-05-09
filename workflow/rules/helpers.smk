# workflow/rules/helpers.smk
import os
import re

SAMPLES = list(config["samples"].keys())
SAMPLE_SET = set(SAMPLES)
OUTDIR = config["output"]["dir"]
KEEP_BAM = bool(config.get("output", {}).get("keep_bam", False))

BWA_INDEX_EXTS = ("amb", "ann", "bwt", "pac", "sa")
BWA_INDEXED_FASTA = config["reference"]["bwa_indexed_fasta"]

CIRI3_NORMALIZE_SCRIPT = os.path.join(
    workflow.basedir,
    "rules",
    "scripts",
    "normalize_ciri3_outputs.py",
)
CIRI3_RATIO_RELATIVE_DEG_SCRIPT = os.path.join(
    workflow.basedir,
    "rules",
    "scripts",
    "compute_ciri3_ratio_relative_deg.R",
)


def _sanitize_path_component(value):
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")
    sanitized = "".join(ch if ch in allowed else "_" for ch in str(value))
    return sanitized.strip("._") or "default"


def _as_bool(value, default=False):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "ture", "1", "yes", "y", "on"}:
            return True
        if normalized in {"false", "0", "no", "n", "off"}:
            return False
    return bool(value)


def _get_bool(cfg, key, default):
    return _as_bool(cfg.get(key, default), default=default)


def _numeric_config(cfg, key, default, cast, label, minimum=None, maximum=None, inclusive_min=True, inclusive_max=True):
    value = cfg.get(key, default)
    try:
        if cast is int:
            if isinstance(value, bool):
                raise ValueError
            value_as_float = float(value)
            if not value_as_float.is_integer():
                raise ValueError
            parsed = int(value_as_float)
        else:
            parsed = cast(value)
    except (TypeError, ValueError):
        raise ValueError(f"{label} must be numeric; got {value!r}.")

    if minimum is not None:
        if inclusive_min and parsed < minimum:
            raise ValueError(f"{label} must be >= {minimum}; got {parsed}.")
        if not inclusive_min and parsed <= minimum:
            raise ValueError(f"{label} must be > {minimum}; got {parsed}.")
    if maximum is not None:
        if inclusive_max and parsed > maximum:
            raise ValueError(f"{label} must be <= {maximum}; got {parsed}.")
        if not inclusive_max and parsed >= maximum:
            raise ValueError(f"{label} must be < {maximum}; got {parsed}.")
    return parsed


ciri3_cfg = config.get("ciri3", {}) or {}
if not isinstance(ciri3_cfg, dict):
    raise ValueError("Top-level ciri3 config must be a mapping.")

CIRI3_REPO_URL = ciri3_cfg.get("repo_url", "https://github.com/gyjames/CIRI3.git")
CIRI3_REPO_REF = str(ciri3_cfg.get("ref", "v3.0.1"))
CIRI3_JAR_NAME = ciri3_cfg.get("jar_name", "CIRI3_Java_18.0.1.jar")
CIRI3_INSTALL_DIR = ciri3_cfg.get(
    "install_dir",
    f"{OUTDIR}/resources/ciri3/{_sanitize_path_component(CIRI3_REPO_REF)}",
)
CIRI3_READY = os.path.join(CIRI3_INSTALL_DIR, ".snakemake_ready")
CIRI3_JAR = os.path.join(CIRI3_INSTALL_DIR, CIRI3_JAR_NAME)
CIRI3_BSJ_YES = os.path.join(CIRI3_INSTALL_DIR, "scripts", "BSJ_yes.R")

rmats_cfg = config.get("rmats_turbo", {}) or {}
if not isinstance(rmats_cfg, dict):
    raise ValueError("Top-level rmats_turbo config must be a mapping.")

RMATS_ENABLED = _get_bool(rmats_cfg, "enabled", False)
RMATS_REPO_URL = rmats_cfg.get("repo_url", "https://github.com/Xinglab/rmats-turbo.git")
RMATS_REPO_REF = str(rmats_cfg.get("ref", "v4.3.0"))
RMATS_INSTALL_DIR = rmats_cfg.get(
    "install_dir",
    f"{OUTDIR}/resources/rmats_turbo/{_sanitize_path_component(RMATS_REPO_REF)}",
)
RMATS_READY = os.path.join(RMATS_INSTALL_DIR, ".snakemake_ready")
RMATS_SCRIPT = os.path.join(RMATS_INSTALL_DIR, "rmats.py")
RMATS_EVENT_TYPES = [str(event_type) for event_type in rmats_cfg.get(
    "event_types",
    ["SE", "A5SS", "A3SS", "MXE", "RI"],
)]
RMATS_VALID_EVENT_TYPES = {"SE", "A5SS", "A3SS", "MXE", "RI"}
unknown_rmats_event_types = sorted(set(RMATS_EVENT_TYPES) - RMATS_VALID_EVENT_TYPES)
if unknown_rmats_event_types:
    raise ValueError(
        "rmats_turbo.event_types contains unsupported event type(s): "
        + ", ".join(unknown_rmats_event_types)
        + ". Valid event types are: A3SS, A5SS, MXE, RI, SE."
    )
RMATS_EVENT_TYPES = [
    event_type
    for event_type in ["SE", "A5SS", "A3SS", "MXE", "RI"]
    if event_type in set(RMATS_EVENT_TYPES)
]
RMATS_EVENT_TYPE_REGEX = "|".join(RMATS_EVENT_TYPES) or r"(?!)"
RMATS_COUNT_TYPES = [str(count_type).upper() for count_type in rmats_cfg.get(
    "count_types",
    ["JC", "JCEC"],
)]
RMATS_VALID_COUNT_TYPES = {"JC", "JCEC"}
unknown_rmats_count_types = sorted(set(RMATS_COUNT_TYPES) - RMATS_VALID_COUNT_TYPES)
if unknown_rmats_count_types:
    raise ValueError(
        "rmats_turbo.count_types contains unsupported count type(s): "
        + ", ".join(unknown_rmats_count_types)
        + ". Valid count types are: JC, JCEC."
    )
RMATS_COUNT_TYPES = [
    count_type for count_type in ["JC", "JCEC"] if count_type in set(RMATS_COUNT_TYPES)
]
RMATS_COUNT_TYPE_REGEX = "|".join(RMATS_COUNT_TYPES) or r"(?!)"
RMATS_READ_TYPE = str(rmats_cfg.get("read_type", "paired"))
if RMATS_READ_TYPE not in {"paired", "single"}:
    raise ValueError("rmats_turbo.read_type must be either 'paired' or 'single'.")
RMATS_LIB_TYPE = str(rmats_cfg.get("lib_type", "fr-unstranded"))
if RMATS_LIB_TYPE not in {"fr-unstranded", "fr-firststrand", "fr-secondstrand"}:
    raise ValueError(
        "rmats_turbo.lib_type must be one of: fr-unstranded, fr-firststrand, fr-secondstrand."
    )
RMATS_READ_LENGTH = _numeric_config(
    rmats_cfg,
    "read_length",
    150,
    int,
    "rmats_turbo.read_length",
    minimum=1,
)
RMATS_TASK = str(rmats_cfg.get("task", "both"))
if RMATS_TASK not in {"prep", "post", "both", "inte", "stat"}:
    raise ValueError("rmats_turbo.task must be one of: prep, post, both, inte, stat.")
RMATS_VARIABLE_READ_LENGTH = _get_bool(rmats_cfg, "variable_read_length", True)
RMATS_ALLOW_CLIPPING = _get_bool(rmats_cfg, "allow_clipping", False)
RMATS_NOVEL_SS = _get_bool(rmats_cfg, "novel_splice_sites", False)
RMATS_PAIRED_STATS = _get_bool(rmats_cfg, "paired_stats", False)
RMATS_STATOFF = _get_bool(rmats_cfg, "statoff", False)
RMATS_INDIVIDUAL_COUNTS = _get_bool(rmats_cfg, "individual_counts", False)
RMATS_EXTRA_ARGS = str(rmats_cfg.get("extra_args", "") or "")

sashimi_cfg = config.get("sashimi", {}) or {}
if not isinstance(sashimi_cfg, dict):
    raise ValueError("Top-level sashimi config must be a mapping.")

SASHIMI_ENABLED = _get_bool(sashimi_cfg, "enabled", False)
SASHIMI_COUNT_TYPE = str(sashimi_cfg.get("count_type", "JC")).upper()
if SASHIMI_COUNT_TYPE not in RMATS_VALID_COUNT_TYPES:
    raise ValueError("sashimi.count_type must be either 'JC' or 'JCEC'.")
if SASHIMI_COUNT_TYPE not in RMATS_COUNT_TYPES:
    RMATS_COUNT_TYPES.append(SASHIMI_COUNT_TYPE)
    RMATS_COUNT_TYPES = [
        count_type for count_type in ["JC", "JCEC"] if count_type in set(RMATS_COUNT_TYPES)
    ]
    RMATS_COUNT_TYPE_REGEX = "|".join(RMATS_COUNT_TYPES) or r"(?!)"
RMATS2SASHIMI_REPO_URL = sashimi_cfg.get(
    "repo_url",
    "https://github.com/Xinglab/rmats2sashimiplot.git",
)
RMATS2SASHIMI_REPO_REF = str(sashimi_cfg.get("ref", "v3.0.0"))
RMATS2SASHIMI_INSTALL_DIR = sashimi_cfg.get(
    "install_dir",
    f"{OUTDIR}/resources/rmats2sashimiplot/{_sanitize_path_component(RMATS2SASHIMI_REPO_REF)}",
)
RMATS2SASHIMI_READY = os.path.join(RMATS2SASHIMI_INSTALL_DIR, ".snakemake_ready")
RMATS2SASHIMI_SCRIPT = os.path.join(
    RMATS2SASHIMI_INSTALL_DIR,
    "src",
    "rmats2sashimiplot",
    "rmats2sashimiplot.py",
)
SASHIMI_FDR_CUTOFF = _numeric_config(
    sashimi_cfg,
    "fdr_cutoff",
    0.05,
    float,
    "sashimi.fdr_cutoff",
    minimum=0,
    maximum=1,
)
SASHIMI_PVALUE_CUTOFF = _numeric_config(
    sashimi_cfg,
    "pvalue_cutoff",
    1.0,
    float,
    "sashimi.pvalue_cutoff",
    minimum=0,
    maximum=1,
)
SASHIMI_INC_DIFF_CUTOFF = _numeric_config(
    sashimi_cfg,
    "inc_diff_cutoff",
    0.1,
    float,
    "sashimi.inc_diff_cutoff",
    minimum=0,
)
SASHIMI_MAX_EVENTS = _numeric_config(
    sashimi_cfg,
    "max_events_per_comparison_type",
    0,
    int,
    "sashimi.max_events_per_comparison_type",
    minimum=0,
)
SASHIMI_MIN_COUNTS = _numeric_config(
    sashimi_cfg,
    "min_counts",
    0,
    int,
    "sashimi.min_counts",
    minimum=0,
)
SASHIMI_EXON_SCALE = _numeric_config(
    sashimi_cfg,
    "exon_s",
    1,
    int,
    "sashimi.exon_s",
    minimum=1,
)
SASHIMI_INTRON_SCALE = _numeric_config(
    sashimi_cfg,
    "intron_s",
    5,
    int,
    "sashimi.intron_s",
    minimum=1,
)
SASHIMI_FONT_SIZE = _numeric_config(
    sashimi_cfg,
    "font_size",
    8,
    int,
    "sashimi.font_size",
    minimum=1,
)
SASHIMI_FIG_WIDTH = _numeric_config(
    sashimi_cfg,
    "fig_width",
    8,
    float,
    "sashimi.fig_width",
    minimum=0,
    inclusive_min=False,
)
SASHIMI_FIG_HEIGHT = _numeric_config(
    sashimi_cfg,
    "fig_height",
    0,
    float,
    "sashimi.fig_height",
    minimum=0,
)
SASHIMI_FAIL_ON_ERROR = _get_bool(sashimi_cfg, "fail_on_error", False)
SASHIMI_KEEP_EVENT_CHR_PREFIX = _get_bool(sashimi_cfg, "keep_event_chr_prefix", False)
SASHIMI_REMOVE_EVENT_CHR_PREFIX = _get_bool(sashimi_cfg, "remove_event_chr_prefix", False)
if SASHIMI_KEEP_EVENT_CHR_PREFIX and SASHIMI_REMOVE_EVENT_CHR_PREFIX:
    raise ValueError(
        "sashimi.keep_event_chr_prefix and sashimi.remove_event_chr_prefix cannot both be true."
    )
SASHIMI_COLORS = sashimi_cfg.get("colors", [])
if isinstance(SASHIMI_COLORS, str):
    SASHIMI_COLORS = [color.strip() for color in SASHIMI_COLORS.split(",") if color.strip()]
else:
    SASHIMI_COLORS = [str(color) for color in SASHIMI_COLORS]
SASHIMI_EXTRA_ARGS = str(sashimi_cfg.get("extra_args", "") or "")


def maybe_temp(path):
    return path if KEEP_BAM else temp(path)


def _auto_comparisons_from_groups(groups):
    group_names = list(groups.keys())
    comparisons = {}
    for idx in range(1, len(group_names)):
        case_group = group_names[idx]
        control_group = group_names[idx - 1]
        comparisons[f"{case_group}_vs_{control_group}"] = {
            "case_group": case_group,
            "control_group": control_group,
            "case": list(groups[case_group]),
            "control": list(groups[control_group]),
        }
    for gap in range(2, len(group_names)):
        for idx in range(gap, len(group_names)):
            case_group = group_names[idx]
            control_group = group_names[idx - gap]
            comparisons[f"{case_group}_vs_{control_group}"] = {
                "case_group": case_group,
                "control_group": control_group,
                "case": list(groups[case_group]),
                "control": list(groups[control_group]),
            }
    return comparisons


def _format_valid_samples():
    return ", ".join(SAMPLES) if SAMPLES else "(none configured)"


def _unknown_sample_error(context, samples):
    unknown = sorted(set(samples) - SAMPLE_SET)
    if unknown:
        raise ValueError(
            f"{context} includes unknown sample names: {', '.join(unknown)}. "
            f"Valid top-level samples are: {_format_valid_samples()}."
        )


def _comparison_from_spec(name, spec, groups):
    if not isinstance(spec, dict):
        raise ValueError(f"deg.comparisons.{name} must be a mapping.")

    def resolve_group_or_samples(role):
        group_key = f"{role}_group"
        value = spec.get(group_key, spec.get(role))
        if value is None:
            raise ValueError(
                f"deg.comparisons.{name} must define '{role}_group' or '{role}'."
            )

        if isinstance(value, str):
            if value not in groups:
                raise ValueError(
                    f"deg.comparisons.{name}.{role}_group references unknown group "
                    f"'{value}'. Valid groups are: {', '.join(groups.keys())}."
                )
            return value, list(groups[value])

        try:
            samples = list(value)
        except TypeError:
            raise ValueError(
                f"deg.comparisons.{name}.{role} must be a group name or a list of samples."
            )
        _unknown_sample_error(f"deg.comparisons.{name}.{role}", samples)
        sample_set = set(samples)
        matches = [
            group_name
            for group_name, members in groups.items()
            if set(members) == sample_set
        ]
        if len(matches) != 1:
            raise ValueError(
                f"deg.comparisons.{name}.{role} must match exactly one deg.groups entry. "
                "Use group names in comparisons, or update deg.groups so it is the single "
                "source of truth for sample membership."
            )
        return matches[0], samples

    case_group, case_samples = resolve_group_or_samples("case")
    control_group, control_samples = resolve_group_or_samples("control")
    if case_group == control_group:
        raise ValueError(f"deg.comparisons.{name} compares group '{case_group}' to itself.")
    return {
        "case_group": case_group,
        "control_group": control_group,
        "case": case_samples,
        "control": control_samples,
    }


def _normalize_explicit_comparisons(comparisons, groups):
    normalized = {}
    for name, spec in comparisons.items():
        if not name or "/" in str(name):
            raise ValueError(
                f"DEG comparison name {name!r} is invalid; names are used in output paths "
                "and cannot be empty or contain '/'."
            )
        normalized[str(name)] = _comparison_from_spec(str(name), spec, groups)
    return normalized


deg_cfg = config.get("deg", {}) or {}
if not isinstance(deg_cfg, dict):
    raise ValueError("Top-level deg config must be a mapping.")

if "enabled" in deg_cfg:
    raise ValueError("Deprecated 'deg.enabled' is no longer supported. Use 'deg.run_deseq2' instead.")

if "ciri3_de" in config:
    raise ValueError(
        "Deprecated top-level 'ciri3_de:' is no longer supported. Move CIRI3 "
        "differential-expression settings into the unified top-level 'deg:' block."
    )

DEG_RUN_DESEQ2 = _get_bool(deg_cfg, "run_deseq2", False)
CIRI3_DE_RUN_BSJ = _get_bool(deg_cfg, "run_de_bsj", False)
CIRI3_DE_RUN_RATIO = _get_bool(deg_cfg, "run_de_ratio", False)
CIRI3_DE_RUN_RELATIVE = _get_bool(deg_cfg, "run_de_relative", False)
CIRI3_DE_ENABLED = CIRI3_DE_RUN_BSJ or CIRI3_DE_RUN_RATIO or CIRI3_DE_RUN_RELATIVE

DEG_GROUPS = deg_cfg.get("groups", {}) or {}
if DEG_GROUPS and not isinstance(DEG_GROUPS, dict):
    raise ValueError("deg.groups must be a group-name-to-sample-list mapping.")
DEG_GROUP_NAMES = list(DEG_GROUPS.keys())

DEG_MIN_TOTAL_COUNT = _numeric_config(
    deg_cfg, "min_total_count", 10, int, "deg.min_total_count", minimum=0
)
DEG_MIN_SAMPLES_DETECTED = _numeric_config(
    deg_cfg, "min_samples_detected", 2, int, "deg.min_samples_detected", minimum=1
)
DEG_PADJ_CUTOFF = _numeric_config(
    deg_cfg, "padj_cutoff", 0.05, float, "deg.padj_cutoff", minimum=0, maximum=1,
    inclusive_min=False, inclusive_max=False
)
DEG_LFC_CUTOFF = _numeric_config(
    deg_cfg, "lfc_cutoff", 1.0, float, "deg.lfc_cutoff", minimum=0
)

CIRI3_DE_OUTDIR = deg_cfg.get("ciri3_outdir", f"{OUTDIR}/deg/ciri3")
CIRI3_DE_USE_FEATURECOUNTS = _get_bool(
    deg_cfg,
    "ciri3_gene_expression_from_featurecounts",
    True,
)

if SASHIMI_ENABLED and not RMATS_ENABLED:
    raise ValueError("sashimi.enabled requires rmats_turbo.enabled because plots use rMATS event files.")

RMATS_WORKFLOW_ACTIVE = RMATS_ENABLED or SASHIMI_ENABLED
DEG_ACTIVE = DEG_RUN_DESEQ2 or CIRI3_DE_ENABLED or RMATS_WORKFLOW_ACTIVE

if DEG_ACTIVE:
    if len(DEG_GROUP_NAMES) < 2:
        raise ValueError("DEG/CIRI3/rMATS analysis is enabled, but fewer than 2 deg.groups are configured.")

    sample_to_group = {}
    duplicated = set()
    for group_name, members in DEG_GROUPS.items():
        if not group_name or "/" in str(group_name):
            raise ValueError(
                f"DEG group name {group_name!r} is invalid; group names are used in output "
                "paths and cannot be empty or contain '/'."
            )
        if isinstance(members, str) or not hasattr(members, "__iter__"):
            raise ValueError(f"deg.groups.{group_name} must be a list of sample names.")
        members = list(members)
        if len(members) < 2:
            raise ValueError(
                f"deg.groups.{group_name} has {len(members)} sample(s). At least two "
                "biological replicates per group are required for enabled comparison modules."
            )
        _unknown_sample_error(f"deg.groups.{group_name}", members)
        for sample in members:
            if sample in sample_to_group:
                duplicated.add(sample)
            sample_to_group[sample] = group_name
    if duplicated:
        raise ValueError(
            "Each sample can only belong to one DEG group. Duplicated samples: "
            + ", ".join(sorted(duplicated))
        )

    selected_deg_samples = [
        sample
        for members in DEG_GROUPS.values()
        for sample in members
    ]
    if DEG_MIN_SAMPLES_DETECTED > len(selected_deg_samples):
        raise ValueError(
            "deg.min_samples_detected cannot exceed the number of samples in deg.groups "
            f"({len(selected_deg_samples)})."
        )

    if CIRI3_DE_RUN_BSJ and not CIRI3_DE_USE_FEATURECOUNTS:
        raise ValueError(
            "deg.run_de_bsj requires deg.ciri3_gene_expression_from_featurecounts: true; "
            "no alternate gene-expression source is currently implemented."
        )

explicit_deg_comparisons = deg_cfg.get("comparisons")

if DEG_ACTIVE:
    if explicit_deg_comparisons is not None:
        DEG_COMPARISONS = _normalize_explicit_comparisons(explicit_deg_comparisons, DEG_GROUPS)
    else:
        DEG_COMPARISONS = _auto_comparisons_from_groups(DEG_GROUPS)
else:
    DEG_COMPARISONS = {}

DEG_COMPARISON_NAMES = list(DEG_COMPARISONS.keys())
CIRI3_DE_COMPARISON_NAMES = DEG_COMPARISON_NAMES
CIRI3_DE_COMPARISON_REGEX = "|".join(
    re.escape(name) for name in CIRI3_DE_COMPARISON_NAMES
) or r"(?!)"
RMATS_COMPARISON_NAMES = DEG_COMPARISON_NAMES if RMATS_ENABLED else []
RMATS_COMPARISON_REGEX = "|".join(
    re.escape(name) for name in RMATS_COMPARISON_NAMES
) or r"(?!)"

SASHIMI_LABEL_B1_CONFIG = sashimi_cfg.get("label_b1")
SASHIMI_LABEL_B2_CONFIG = sashimi_cfg.get("label_b2")


def star_bam_path(sample):
    return f"{OUTDIR}/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam"


def star_bai_path(sample):
    return f"{OUTDIR}/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"


def rmats_comparison_case_samples(comparison):
    return list(DEG_COMPARISONS[comparison]["case"])


def rmats_comparison_control_samples(comparison):
    return list(DEG_COMPARISONS[comparison]["control"])


def rmats_sample_bams(samples):
    return [star_bam_path(sample) for sample in samples]


def rmats_sample_bais(samples):
    return [star_bai_path(sample) for sample in samples]


def rmats_comparison_bams(wildcards):
    samples = rmats_comparison_case_samples(wildcards.comparison)
    samples += rmats_comparison_control_samples(wildcards.comparison)
    return rmats_sample_bams(samples)


def rmats_comparison_bais(wildcards):
    samples = rmats_comparison_case_samples(wildcards.comparison)
    samples += rmats_comparison_control_samples(wildcards.comparison)
    return rmats_sample_bais(samples)


def sashimi_comparison_group_names(comparison):
    comp = DEG_COMPARISONS[comparison]
    return [comp["case_group"], comp["control_group"]]


def sashimi_comparison_groups(comparison):
    comp = DEG_COMPARISONS[comparison]
    return {
        comp["case_group"]: list(comp["case"]),
        comp["control_group"]: list(comp["control"]),
    }


def _configured_or_default_label(configured, default):
    return str(configured) if configured not in (None, "") else str(default)


def sashimi_label_b1(wildcards):
    return _configured_or_default_label(
        SASHIMI_LABEL_B1_CONFIG,
        DEG_COMPARISONS[wildcards.comparison]["case_group"],
    )


def sashimi_label_b2(wildcards):
    return _configured_or_default_label(
        SASHIMI_LABEL_B2_CONFIG,
        DEG_COMPARISONS[wildcards.comparison]["control_group"],
    )


RMATS_EVENT_TARGETS = (
    [
        f"{OUTDIR}/rmats/{comparison}/{event_type}.MATS.{count_type}.txt"
        for comparison in RMATS_COMPARISON_NAMES
        for event_type in RMATS_EVENT_TYPES
        for count_type in RMATS_COUNT_TYPES
    ]
    if RMATS_ENABLED
    else []
)
SASHIMI_TARGETS = (
    [
        f"{OUTDIR}/rmats/sashimi/{comparison}/{event_type}.{SASHIMI_COUNT_TYPE}/manifest.tsv"
        for comparison in RMATS_COMPARISON_NAMES
        for event_type in RMATS_EVENT_TYPES
    ]
    + [
        f"{OUTDIR}/rmats/sashimi/{comparison}/{event_type}.{SASHIMI_COUNT_TYPE}/plots.done"
        for comparison in RMATS_COMPARISON_NAMES
        for event_type in RMATS_EVENT_TYPES
    ]
    if SASHIMI_ENABLED
    else []
)
CIRI3_DE_METHODS = ["de_bsj", "de_ratio", "de_relative"]
CIRI3_DE_ENABLED_METHODS = (
    (["de_bsj"] if CIRI3_DE_RUN_BSJ else [])
    + (["de_ratio"] if CIRI3_DE_RUN_RATIO else [])
    + (["de_relative"] if CIRI3_DE_RUN_RELATIVE else [])
)

motif_cfg = config.get("motif", {})
MOTIF_ENABLED = _as_bool(motif_cfg.get("enabled", True), default=True)

if DEG_ACTIVE and not DEG_COMPARISON_NAMES:
    raise ValueError("DEG/CIRI3/rMATS analysis is enabled, but no comparisons could be generated.")

if DEG_ACTIVE:
    for comp_name, comp in DEG_COMPARISONS.items():
        case_samples = list(comp.get("case", []))
        control_samples = list(comp.get("control", []))
        if len(case_samples) < 2 or len(control_samples) < 2:
            raise ValueError(
                f"DEG comparison '{comp_name}' must have at least two biological "
                "replicates in both case and control groups."
            )
        overlap = sorted(set(case_samples) & set(control_samples))
        if overlap:
            raise ValueError(
                f"DEG comparison '{comp_name}' has samples in both case and control: "
                + ", ".join(overlap)
            )


def ciri3_de_case_samples(comparison):
    return list(DEG_COMPARISONS[comparison]["case"])


def ciri3_de_control_samples(comparison):
    return list(DEG_COMPARISONS[comparison]["control"])


def ciri3_de_all_samples(comparison):
    return ciri3_de_case_samples(comparison) + ciri3_de_control_samples(comparison)


def ciri3_de_sample_class_map(comparison):
    out = {}
    for sample in ciri3_de_case_samples(comparison):
        out[sample] = "Case"
    for sample in ciri3_de_control_samples(comparison):
        out[sample] = "Control"
    return out


def ciri3_per_sample_result_paths(comparison):
    return [
        f"{OUTDIR}/ciri3/per_sample/{sample}.ciri3"
        for sample in ciri3_de_all_samples(comparison)
    ]
