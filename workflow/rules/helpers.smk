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

# CIRI3 arm master switch. Defaults to True so existing configs that pre-date
# circtools integration keep running the CIRI3 detection + downstream BSJ/FSJ
# pipeline unchanged. Set ciri3.enabled: false to skip CIRI3 entirely when
# running an exclusively-circtools-driven workflow (which then requires
# circtools.detect.enabled to provide a canonical DCC).
CIRI3_ENABLED = _get_bool(ciri3_cfg, "enabled", True)

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
CIRI3_RMATS_EXE = os.path.join(CIRI3_INSTALL_DIR, "scripts", "rMATSexe")

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
RMATS2SASHIMI_READY = os.path.join(RMATS2SASHIMI_INSTALL_DIR, ".snakemake_py3_miso_ready")
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
SASHIMI_AUTO_SCALE = _get_bool(sashimi_cfg, "auto_scale", True)
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
SASHIMI_MIN_FIG_WIDTH = _numeric_config(
    sashimi_cfg,
    "min_fig_width",
    SASHIMI_FIG_WIDTH,
    float,
    "sashimi.min_fig_width",
    minimum=0,
    inclusive_min=False,
)
SASHIMI_MAX_FIG_WIDTH = _numeric_config(
    sashimi_cfg,
    "max_fig_width",
    max(18, SASHIMI_MIN_FIG_WIDTH),
    float,
    "sashimi.max_fig_width",
    minimum=SASHIMI_MIN_FIG_WIDTH,
)
SASHIMI_MIN_FIG_HEIGHT = _numeric_config(
    sashimi_cfg,
    "min_fig_height",
    SASHIMI_FIG_HEIGHT if SASHIMI_FIG_HEIGHT > 0 else 4.5,
    float,
    "sashimi.min_fig_height",
    minimum=0,
    inclusive_min=False,
)
SASHIMI_MAX_FIG_HEIGHT = _numeric_config(
    sashimi_cfg,
    "max_fig_height",
    max(12, SASHIMI_MIN_FIG_HEIGHT),
    float,
    "sashimi.max_fig_height",
    minimum=SASHIMI_MIN_FIG_HEIGHT,
)
SASHIMI_MIN_FONT_SIZE = _numeric_config(
    sashimi_cfg,
    "min_font_size",
    min(6, SASHIMI_FONT_SIZE),
    int,
    "sashimi.min_font_size",
    minimum=1,
    maximum=SASHIMI_FONT_SIZE,
)
SASHIMI_BSJ_FLANK = _numeric_config(
    sashimi_cfg,
    "bsj_flank",
    250,
    int,
    "sashimi.bsj_flank",
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
SASHIMI_GFF3_CONFIG = (
    config.get("reference", {}).get("gff3")
    or sashimi_cfg.get("annotation_gff3")
    or sashimi_cfg.get("gff3")
)
SASHIMI_GFF3_FROM_GTF = SASHIMI_GFF3_CONFIG in (None, "")
SASHIMI_GFF3 = (
    f"{OUTDIR}/resources/annotation/"
    f"{_sanitize_path_component(os.path.basename(config['reference']['gtf']))}.gff3"
    if SASHIMI_GFF3_FROM_GTF
    else str(SASHIMI_GFF3_CONFIG)
)


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

RMATS_WORKFLOW_ACTIVE = RMATS_ENABLED or SASHIMI_ENABLED
DEG_ACTIVE = DEG_RUN_DESEQ2 or CIRI3_DE_ENABLED or RMATS_WORKFLOW_ACTIVE

if SASHIMI_ENABLED and not (DEG_RUN_DESEQ2 or CIRI3_DE_ENABLED):
    raise ValueError(
        "sashimi.enabled requires at least one BSJ differential source: "
        "deg.run_deseq2, deg.run_de_bsj, deg.run_de_ratio, or deg.run_de_relative."
    )

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

    DEG_SELECTED_SAMPLES = [
        sample
        for members in DEG_GROUPS.values()
        for sample in members
    ]
    if DEG_MIN_SAMPLES_DETECTED > len(DEG_SELECTED_SAMPLES):
        raise ValueError(
            "deg.min_samples_detected cannot exceed the number of samples in deg.groups "
            f"({len(DEG_SELECTED_SAMPLES)})."
        )

    if CIRI3_DE_RUN_BSJ and not CIRI3_DE_USE_FEATURECOUNTS:
        raise ValueError(
            "deg.run_de_bsj requires deg.ciri3_gene_expression_from_featurecounts: true; "
            "no alternate gene-expression source is currently implemented."
        )
else:
    DEG_SELECTED_SAMPLES = []

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
RMATS_GROUP_NAMES = DEG_GROUP_NAMES if RMATS_ENABLED else []
RMATS_GROUP_REGEX = "|".join(
    re.escape(name) for name in RMATS_GROUP_NAMES
) or r"(?!)"
RMATS_COMPARISON_NAMES = DEG_COMPARISON_NAMES if RMATS_ENABLED else []
RMATS_COMPARISON_REGEX = "|".join(
    re.escape(name) for name in RMATS_COMPARISON_NAMES
) or r"(?!)"

BSJ_SASHIMI_METHODS = (
    (["deseq2"] if DEG_RUN_DESEQ2 else [])
    + (["de_bsj"] if CIRI3_DE_RUN_BSJ else [])
    + (["de_ratio"] if CIRI3_DE_RUN_RATIO else [])
    + (["de_relative"] if CIRI3_DE_RUN_RELATIVE else [])
)
BSJ_SASHIMI_METHOD_REGEX = "|".join(
    re.escape(method) for method in BSJ_SASHIMI_METHODS
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


def rmats_group_samples(group):
    return list(DEG_GROUPS[group])


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


def rmats_group_bams(wildcards):
    return rmats_sample_bams(rmats_group_samples(wildcards.group))


def rmats_group_bais(wildcards):
    return rmats_sample_bais(rmats_group_samples(wildcards.group))


def sashimi_all_group_bams(wildcards):
    return rmats_sample_bams(DEG_SELECTED_SAMPLES)


def sashimi_all_group_bais(wildcards):
    return rmats_sample_bais(DEG_SELECTED_SAMPLES)


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


def bsj_sashimi_result_path(wildcards):
    if wildcards.method == "deseq2":
        return f"{OUTDIR}/deg/deseq2/pairwise/{wildcards.comparison}/deseq2_results.tsv"
    if wildcards.method in {"de_bsj", "de_ratio", "de_relative"}:
        return f"{CIRI3_DE_OUTDIR}/{wildcards.comparison}/{wildcards.method}/result.txt"
    raise ValueError(f"Unsupported BSJ sashimi method: {wildcards.method}")


RMATS_EVENT_TARGETS = (
    [
        f"{OUTDIR}/rmats/groups/{group}/{event_type}.MATS.{count_type}.txt"
        for group in RMATS_GROUP_NAMES
        for event_type in RMATS_EVENT_TYPES
        for count_type in RMATS_COUNT_TYPES
    ]
    if RMATS_ENABLED
    else []
)
SASHIMI_TARGETS = (
    [
        f"{OUTDIR}/rmats/sashimi/bsj/{method}/{comparison}/manifest.tsv"
        for method in BSJ_SASHIMI_METHODS
        for comparison in DEG_COMPARISON_NAMES
    ]
    + [
        f"{OUTDIR}/rmats/sashimi/bsj/{method}/{comparison}/plots.done"
        for method in BSJ_SASHIMI_METHODS
        for comparison in DEG_COMPARISON_NAMES
    ]
    + [
        f"{OUTDIR}/rmats/sashimi/bsj/{method}/{comparison}/plots_bsj_only"
        for method in BSJ_SASHIMI_METHODS
        for comparison in DEG_COMPARISON_NAMES
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

# ---------------------------------------------------------------------------
# circrna_plots: rich CIRI3 DE_BSJ plots + dataset-wide BSJ summaries.
# ---------------------------------------------------------------------------
circrna_plots_cfg = config.get("circrna_plots", {}) or {}
if not isinstance(circrna_plots_cfg, dict):
    raise ValueError("Top-level circrna_plots config must be a mapping.")

CIRCRNA_PLOTS_ENABLED = _get_bool(circrna_plots_cfg, "enabled", False)

if CIRCRNA_PLOTS_ENABLED:
    if not CIRI3_DE_RUN_BSJ:
        raise ValueError(
            "circrna_plots.enabled requires deg.run_de_bsj: true (CIRI3 DE_BSJ "
            "must produce result.txt + infor.tsv per comparison)."
        )

    CIRCRNA_PLOTS_REVERSE_LFC = _get_bool(
        circrna_plots_cfg, "reverse_ciri3_logfc_direction", True
    )
    CIRCRNA_PLOTS_LABEL_TOP_N = _numeric_config(
        circrna_plots_cfg, "label_top_n", 15, int,
        "circrna_plots.label_top_n", minimum=0,
    )
    CIRCRNA_PLOTS_FDR_FLOOR = _numeric_config(
        circrna_plots_cfg, "fdr_floor_for_plot", 1e-300, float,
        "circrna_plots.fdr_floor_for_plot", minimum=0, inclusive_min=False,
    )
    CIRCRNA_PLOTS_PER_MILLION = _numeric_config(
        circrna_plots_cfg, "per_million", 1e6, float,
        "circrna_plots.per_million", minimum=0, inclusive_min=False,
    )

    raw_top_n = circrna_plots_cfg.get("heatmap_top_n_values", [50, 25])
    if not isinstance(raw_top_n, (list, tuple)) or len(raw_top_n) == 0:
        raise ValueError("circrna_plots.heatmap_top_n_values must be a non-empty list.")
    parsed_top_n = []
    for value in raw_top_n:
        if isinstance(value, bool):
            raise ValueError("circrna_plots.heatmap_top_n_values must be a list of positive integers.")
        try:
            n = int(value)
        except (TypeError, ValueError):
            raise ValueError("circrna_plots.heatmap_top_n_values must be a list of positive integers.")
        if n <= 0:
            raise ValueError("circrna_plots.heatmap_top_n_values entries must be > 0.")
        parsed_top_n.append(n)
    if len(set(parsed_top_n)) != len(parsed_top_n):
        raise ValueError("circrna_plots.heatmap_top_n_values must not contain duplicates.")
    CIRCRNA_PLOTS_HEATMAP_TOP_N_VALUES = parsed_top_n

    CIRCRNA_PLOTS_HEATMAP_CLUSTER_ROWS = _get_bool(circrna_plots_cfg, "heatmap_cluster_rows", True)
    CIRCRNA_PLOTS_HEATMAP_CLUSTER_COLS = _get_bool(circrna_plots_cfg, "heatmap_cluster_cols", False)
    CIRCRNA_PLOTS_HEATMAP_SHOW_ROWNAMES_ALL = _get_bool(
        circrna_plots_cfg, "heatmap_show_rownames_all", False
    )
    CIRCRNA_PLOTS_HEATMAP_SHOW_ROWNAMES_TOP = _get_bool(
        circrna_plots_cfg, "heatmap_show_rownames_top", True
    )
    CIRCRNA_PLOTS_HEATMAP_USE_LOG_COLORMAP = _get_bool(
        circrna_plots_cfg, "heatmap_absolute_use_log_for_colormap", False
    )
    CIRCRNA_PLOTS_HEATMAP_CLIP_Q = _numeric_config(
        circrna_plots_cfg, "heatmap_absolute_upper_clip_quantile", 0.99, float,
        "circrna_plots.heatmap_absolute_upper_clip_quantile",
        minimum=0, maximum=1, inclusive_min=False,
    )
    CIRCRNA_PLOTS_HEATMAP_LOG_PSEUDOCOUNT = _numeric_config(
        circrna_plots_cfg, "heatmap_log_pseudocount", 1, float,
        "circrna_plots.heatmap_log_pseudocount", minimum=0,
    )
    CIRCRNA_PLOTS_HEATMAP_ZSCORE_BREAKS = _numeric_config(
        circrna_plots_cfg, "heatmap_zscore_breaks_limit", 3, float,
        "circrna_plots.heatmap_zscore_breaks_limit", minimum=0, inclusive_min=False,
    )

    # Condition order is the order of deg.groups keys (preserved insertion order).
    CIRCRNA_PLOTS_CONDITION_ORDER = list(DEG_GROUP_NAMES)

    # Auto-pick colors from a Tableau-10-like default palette in DEG_GROUPS order,
    # then allow user override via circrna_plots.condition_colors.
    _CIRCRNA_PLOTS_DEFAULT_PALETTE = [
        "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
        "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
    ]
    _default_colors = {
        name: _CIRCRNA_PLOTS_DEFAULT_PALETTE[i % len(_CIRCRNA_PLOTS_DEFAULT_PALETTE)]
        for i, name in enumerate(CIRCRNA_PLOTS_CONDITION_ORDER)
    }
    _color_overrides = circrna_plots_cfg.get("condition_colors") or {}
    if not isinstance(_color_overrides, dict):
        raise ValueError(
            "circrna_plots.condition_colors must be a mapping of group name -> hex color."
        )
    unknown_color_groups = sorted(set(_color_overrides) - set(CIRCRNA_PLOTS_CONDITION_ORDER))
    if unknown_color_groups:
        raise ValueError(
            "circrna_plots.condition_colors references unknown deg.groups: "
            + ", ".join(unknown_color_groups)
        )
    CIRCRNA_PLOTS_CONDITION_COLORS_RESOLVED = {
        name: str(_color_overrides.get(name, _default_colors[name]))
        for name in CIRCRNA_PLOTS_CONDITION_ORDER
    }

    CIRCRNA_PLOTS_SAMPLES = list(DEG_SELECTED_SAMPLES)
    CIRCRNA_PLOTS_OUTDIR = f"{OUTDIR}/circrna_plots"
else:
    CIRCRNA_PLOTS_REVERSE_LFC = False
    CIRCRNA_PLOTS_LABEL_TOP_N = 15
    CIRCRNA_PLOTS_FDR_FLOOR = 1e-300
    CIRCRNA_PLOTS_PER_MILLION = 1e6
    CIRCRNA_PLOTS_HEATMAP_TOP_N_VALUES = []
    CIRCRNA_PLOTS_HEATMAP_CLUSTER_ROWS = True
    CIRCRNA_PLOTS_HEATMAP_CLUSTER_COLS = False
    CIRCRNA_PLOTS_HEATMAP_SHOW_ROWNAMES_ALL = False
    CIRCRNA_PLOTS_HEATMAP_SHOW_ROWNAMES_TOP = True
    CIRCRNA_PLOTS_HEATMAP_USE_LOG_COLORMAP = False
    CIRCRNA_PLOTS_HEATMAP_CLIP_Q = 0.99
    CIRCRNA_PLOTS_HEATMAP_LOG_PSEUDOCOUNT = 1
    CIRCRNA_PLOTS_HEATMAP_ZSCORE_BREAKS = 3
    CIRCRNA_PLOTS_CONDITION_ORDER = []
    CIRCRNA_PLOTS_CONDITION_COLORS_RESOLVED = {}
    CIRCRNA_PLOTS_SAMPLES = []
    CIRCRNA_PLOTS_OUTDIR = f"{OUTDIR}/circrna_plots"

# ---------------------------------------------------------------------------
# circtools integration (https://github.com/jakobilab/circtools).
#
# All circtools state is parsed here so circtools.smk only references already-
# resolved constants and helper functions. The arm has four sub-modules
# (detect / circtest / primex / quickcheck), each gated by its own enable
# switch under the top-level `circtools:` config block; comparisons follow
# the same explicit-or-auto-from-groups pattern used by `deg.comparisons`.
# ---------------------------------------------------------------------------
circtools_cfg = config.get("circtools", {}) or {}
if not isinstance(circtools_cfg, dict):
    raise ValueError("Top-level circtools config must be a mapping.")

CIRCTOOLS_ENABLED = _get_bool(circtools_cfg, "enabled", False)


def _circtools_sub_cfg(key):
    sub = circtools_cfg.get(key, {}) or {}
    if not isinstance(sub, dict):
        raise ValueError(f"circtools.{key} must be a mapping.")
    return sub


_circtools_detect_cfg = _circtools_sub_cfg("detect")
_circtools_circtest_cfg = _circtools_sub_cfg("circtest")
_circtools_primex_cfg = _circtools_sub_cfg("primex")
_circtools_quickcheck_cfg = _circtools_sub_cfg("quickcheck")
_circtools_star_cfg = _circtools_sub_cfg("star")

CIRCTOOLS_DETECT_ENABLED = (
    CIRCTOOLS_ENABLED and _get_bool(_circtools_detect_cfg, "enabled", False)
)
CIRCTOOLS_CIRCTEST_ENABLED = (
    CIRCTOOLS_ENABLED and _get_bool(_circtools_circtest_cfg, "enabled", False)
)
CIRCTOOLS_PRIMEX_ENABLED = (
    CIRCTOOLS_ENABLED and _get_bool(_circtools_primex_cfg, "enabled", False)
)
CIRCTOOLS_QUICKCHECK_ENABLED = (
    CIRCTOOLS_ENABLED and _get_bool(_circtools_quickcheck_cfg, "enabled", False)
)

# Convenience: any "downstream" circtools sub-module (i.e. needing a canonical
# DCC). detect itself is what PRODUCES the DCC in the detect path so it's not
# a downstream consumer; the three below are.
_CIRCTOOLS_DOWNSTREAM_ENABLED = (
    CIRCTOOLS_CIRCTEST_ENABLED
    or CIRCTOOLS_PRIMEX_ENABLED
    or CIRCTOOLS_QUICKCHECK_ENABLED
)

# Sub-module sanity: if the arm is on, at least one sub-module must be on too.
if CIRCTOOLS_ENABLED and not (
    CIRCTOOLS_DETECT_ENABLED or _CIRCTOOLS_DOWNSTREAM_ENABLED
):
    raise ValueError(
        "circtools.enabled is true but no sub-module is enabled. Set at least "
        "one of circtools.detect.enabled / circtools.circtest.enabled / "
        "circtools.primex.enabled / circtools.quickcheck.enabled."
    )

# Canonical DCC source: any downstream sub-module needs ONE of detect or
# ciri3 to produce a canonical DCC. The detect arm produces a real DCC;
# without it, we synthesize a DCC from CIRI3 BSJ/FSJ matrices, which means
# ciri3 must be enabled.
if _CIRCTOOLS_DOWNSTREAM_ENABLED and not (CIRCTOOLS_DETECT_ENABLED or CIRI3_ENABLED):
    raise ValueError(
        "circtools downstream sub-modules (circtest/primex/quickcheck) require "
        "either circtools.detect.enabled: true OR ciri3.enabled: true to "
        "provide a canonical DCC."
    )

CIRCTOOLS_INSTALL_MARKER = f"{OUTDIR}/resources/circtools/.install_ready"

# Samples used by circtools. Defaults to DEG_SELECTED_SAMPLES (those configured
# in deg.groups) so the per-comparison logic always has a well-defined column
# universe; falls back to all top-level samples when no DEG groups are set.
# This is deliberately a subset of the master SAMPLES list -- circtools detect
# / DCC synthesis runs only over these.
if DEG_SELECTED_SAMPLES:
    CIRCTOOLS_SAMPLES = list(DEG_SELECTED_SAMPLES)
else:
    CIRCTOOLS_SAMPLES = list(SAMPLES)

# ---------------------------------------------------------------------------
# STAR knobs (only used when detect arm is on). Defaults match circtools'
# detect documentation and happen to match the existing CIRI3 STAR rule by
# coincidence, but we keep them independent so circtools' STAR-rule evolution
# does not leak into the CIRI3 path.
# ---------------------------------------------------------------------------
CIRCTOOLS_STAR_CHIM_SEGMENT_MIN = _numeric_config(
    _circtools_star_cfg, "chim_segment_min", 15, int,
    "circtools.star.chim_segment_min", minimum=1,
)
CIRCTOOLS_STAR_CHIM_SCORE_MIN = _numeric_config(
    _circtools_star_cfg, "chim_score_min", 15, int,
    "circtools.star.chim_score_min", minimum=0,
)
CIRCTOOLS_STAR_CHIM_JUNCTION_OVERHANG_MIN = _numeric_config(
    _circtools_star_cfg, "chim_junction_overhang_min", 15, int,
    "circtools.star.chim_junction_overhang_min", minimum=1,
)
CIRCTOOLS_STAR_ALIGN_SJ_OVERHANG_MIN = _numeric_config(
    _circtools_star_cfg, "align_sj_overhang_min", 15, int,
    "circtools.star.align_sj_overhang_min", minimum=1,
)
CIRCTOOLS_STAR_OUT_FILTER_MULTIMAP_NMAX = _numeric_config(
    _circtools_star_cfg, "out_filter_multimap_nmax", 20, int,
    "circtools.star.out_filter_multimap_nmax", minimum=1,
)
CIRCTOOLS_STAR_OUT_FILTER_MISMATCH_NMAX = _numeric_config(
    _circtools_star_cfg, "out_filter_mismatch_nmax", 2, int,
    "circtools.star.out_filter_mismatch_nmax", minimum=0,
)

# ---------------------------------------------------------------------------
# Detect knobs.
# ---------------------------------------------------------------------------
CIRCTOOLS_DETECT_STRANDED = _get_bool(_circtools_detect_cfg, "stranded", True)
CIRCTOOLS_DETECT_SS = _get_bool(_circtools_detect_cfg, "second_strand", False)
CIRCTOOLS_DETECT_FILTER = _get_bool(_circtools_detect_cfg, "filter", True)
# -M: remove circRNA candidates on the mitochondrial chromosome.
CIRCTOOLS_DETECT_CHRM = _get_bool(_circtools_detect_cfg, "remove_chrM", True)
# -fg: also filter by gene annotation (a candidate may not span >1 gene).
CIRCTOOLS_DETECT_FILTER_BY_GENE = _get_bool(
    _circtools_detect_cfg, "filter_by_gene", True
)
# -G + -A + -B: count host gene expression from the sorted BAMs. This is what
# produces a biologically meaningful LinearCount (otherwise LinearCount is
# derived only from junction-spanning linear reads). Strongly recommended.
CIRCTOOLS_DETECT_HOST_GENE = _get_bool(_circtools_detect_cfg, "host_gene", True)
# -R: optional repeats GTF for filtering circRNAs in repetitive regions.
# Empty string (default) means the flag is omitted.
CIRCTOOLS_DETECT_REPEATS_GTF = str(
    _circtools_detect_cfg.get("repeats_gtf", "") or ""
)
if CIRCTOOLS_DETECT_REPEATS_GTF and not os.path.isfile(CIRCTOOLS_DETECT_REPEATS_GTF):
    raise ValueError(
        f"circtools.detect.repeats_gtf is set to '{CIRCTOOLS_DETECT_REPEATS_GTF}' "
        "but that file does not exist. Provide a valid GTF path or remove the "
        "key to skip repeat-region filtering."
    )
CIRCTOOLS_DETECT_COUNT_THRESHOLD = _numeric_config(
    _circtools_detect_cfg, "count_threshold", 2, int,
    "circtools.detect.count_threshold", minimum=0,
)
# Replicate threshold default is 3 (not DCC's library default of 5): with the
# typical 3-replicates-per-group design, a threshold of 5 can never be met and
# would silently filter out every circRNA. Override in config if you have more.
CIRCTOOLS_DETECT_REPLICATE_THRESHOLD = _numeric_config(
    _circtools_detect_cfg, "replicate_threshold", 3, int,
    "circtools.detect.replicate_threshold", minimum=0,
)

# ---------------------------------------------------------------------------
# Circtest knobs.
# ---------------------------------------------------------------------------
CIRCTOOLS_CIRCTEST_MAX_FDR = _numeric_config(
    _circtools_circtest_cfg, "max_fdr", 0.05, float,
    "circtools.circtest.max_fdr",
    minimum=0, maximum=1, inclusive_min=False, inclusive_max=False,
)
CIRCTOOLS_CIRCTEST_PERCENTAGE = _numeric_config(
    _circtools_circtest_cfg, "percentage", 0.01, float,
    "circtools.circtest.percentage",
    minimum=0, maximum=1, inclusive_min=False, inclusive_max=True,
)
CIRCTOOLS_CIRCTEST_FILTER_SAMPLE = _numeric_config(
    _circtools_circtest_cfg, "filter_sample", 3, int,
    "circtools.circtest.filter_sample", minimum=1,
)
CIRCTOOLS_CIRCTEST_FILTER_COUNT = _numeric_config(
    _circtools_circtest_cfg, "filter_count", 5, int,
    "circtools.circtest.filter_count", minimum=0,
)
CIRCTOOLS_CIRCTEST_MAX_PLOTS = _numeric_config(
    _circtools_circtest_cfg, "max_plots", 50, int,
    "circtools.circtest.max_plots", minimum=0,
)

# Circtest comparisons. Explicit comparisons under circtools.circtest.comparisons
# take precedence; otherwise we fall back to DEG_COMPARISONS (which itself was
# derived from deg.comparisons or auto-generated from deg.groups). This keeps
# the configured comparison universe consistent between the DEG arm and
# circtools when the user wants them to match without duplicating config.
_explicit_circtest_comparisons = _circtools_circtest_cfg.get("comparisons")

if _CIRCTOOLS_DOWNSTREAM_ENABLED:
    if _explicit_circtest_comparisons is not None:
        if not DEG_GROUPS:
            raise ValueError(
                "circtools.circtest.comparisons referenced group names but "
                "deg.groups is empty; populate deg.groups (which defines the "
                "valid group->sample mapping for both DEG and circtools)."
            )
        CIRCTOOLS_CIRCTEST_COMPARISONS = _normalize_explicit_comparisons(
            _explicit_circtest_comparisons, DEG_GROUPS
        )
    elif DEG_COMPARISONS:
        # Inherit from the DEG arm (this is the common case).
        CIRCTOOLS_CIRCTEST_COMPARISONS = dict(DEG_COMPARISONS)
    elif DEG_GROUPS:
        # DEG arm wasn't activated but groups are still defined; fall back to
        # the same auto-pairwise convention used elsewhere.
        CIRCTOOLS_CIRCTEST_COMPARISONS = _auto_comparisons_from_groups(DEG_GROUPS)
    else:
        raise ValueError(
            "circtools downstream sub-modules need at least one comparison. "
            "Set either circtools.circtest.comparisons or deg.groups."
        )

    # Validate sample membership: every case+control sample of every
    # comparison must be in CIRCTOOLS_SAMPLES, otherwise the canonical DCC
    # won't contain its column and the subset rule will fail late.
    _circtools_sample_set = set(CIRCTOOLS_SAMPLES)
    for _name, _comp in CIRCTOOLS_CIRCTEST_COMPARISONS.items():
        _missing = sorted(
            (set(_comp["case"]) | set(_comp["control"])) - _circtools_sample_set
        )
        if _missing:
            raise ValueError(
                f"circtools comparison '{_name}' references sample(s) not in "
                f"CIRCTOOLS_SAMPLES: {', '.join(_missing)}. Ensure they appear "
                "in deg.groups (or the explicit comparison's case/control list) "
                "AND in the top-level samples mapping."
            )
        if len(_comp["case"]) < 2 or len(_comp["control"]) < 2:
            raise ValueError(
                f"circtools comparison '{_name}' must have at least 2 "
                "replicates in both case and control; circtest's beta-binomial "
                "ratio test is unstable with fewer."
            )
else:
    CIRCTOOLS_CIRCTEST_COMPARISONS = {}

CIRCTOOLS_CIRCTEST_COMPARISON_NAMES = list(CIRCTOOLS_CIRCTEST_COMPARISONS.keys())
CIRCTOOLS_CIRCTEST_COMPARISON_REGEX = "|".join(
    re.escape(name) for name in CIRCTOOLS_CIRCTEST_COMPARISON_NAMES
) or r"(?!)"

# ---------------------------------------------------------------------------
# Primex knobs.
# ---------------------------------------------------------------------------
CIRCTOOLS_PRIMEX_TOP_N = _numeric_config(
    _circtools_primex_cfg, "top_n", 25, int,
    "circtools.primex.top_n", minimum=1,
)
_primex_size_raw = _circtools_primex_cfg.get("product_size", [80, 160])
if (not isinstance(_primex_size_raw, (list, tuple))
        or len(_primex_size_raw) != 2):
    raise ValueError(
        "circtools.primex.product_size must be a 2-element list [low, high]."
    )
try:
    _primex_size_low = int(_primex_size_raw[0])
    _primex_size_high = int(_primex_size_raw[1])
except (TypeError, ValueError):
    raise ValueError("circtools.primex.product_size values must be integers.")
if _primex_size_low <= 0 or _primex_size_high <= _primex_size_low:
    raise ValueError(
        "circtools.primex.product_size must satisfy 0 < low < high."
    )
CIRCTOOLS_PRIMEX_PRODUCT_SIZE = [_primex_size_low, _primex_size_high]
CIRCTOOLS_PRIMEX_NUM_PAIRS = _numeric_config(
    _circtools_primex_cfg, "num_pairs", 5, int,
    "circtools.primex.num_pairs", minimum=1,
)
CIRCTOOLS_PRIMEX_JUNCTION = str(_circtools_primex_cfg.get("junction", "n"))
if CIRCTOOLS_PRIMEX_JUNCTION not in {"r", "n", "f"}:
    raise ValueError(
        "circtools.primex.junction must be one of: r (reverse primer on BSJ), "
        "n (neither primer on BSJ; default), f (forward primer on BSJ)."
    )
CIRCTOOLS_PRIMEX_NO_BLAST = _get_bool(_circtools_primex_cfg, "no_blast", True)
CIRCTOOLS_PRIMEX_ORGANISM = str(_circtools_primex_cfg.get("organism", "") or "")
if CIRCTOOLS_PRIMEX_ORGANISM and CIRCTOOLS_PRIMEX_ORGANISM not in {"mm", "rn", "hs", "ss"}:
    raise ValueError(
        "circtools.primex.organism must be one of: mm (Mus musculus), "
        "rn (Rattus norvegicus), hs (Homo sapiens), ss (Sus scrofa), or "
        "empty (no BLAST; required for organisms primex doesn't bundle, "
        "including Aedes aegypti)."
    )
CIRCTOOLS_PRIMEX_TITLE = str(_circtools_primex_cfg.get("title", "circtools_primex"))

# ---------------------------------------------------------------------------
# Canonical DCC source. Resolved from arm flags so circtools.smk rules can
# depend on a single set of paths whether they came from real detect or
# CIRI3-synthesized DCC.
# ---------------------------------------------------------------------------
if CIRCTOOLS_DETECT_ENABLED:
    CIRCTOOLS_CANONICAL_DCC_DIR = f"{OUTDIR}/circtools/detect"
else:
    # Implicitly the ciri3-synthesis path (validated above to require CIRI3
    # when no detect and any downstream is on; an unused canonical path is
    # harmless when no downstream is enabled).
    CIRCTOOLS_CANONICAL_DCC_DIR = f"{OUTDIR}/circtools/dcc_from_ciri3/dataset"

CIRCTOOLS_CANONICAL_DCC_CIRC_COUNT = f"{CIRCTOOLS_CANONICAL_DCC_DIR}/CircRNACount"
CIRCTOOLS_CANONICAL_DCC_CIRC_COORD = f"{CIRCTOOLS_CANONICAL_DCC_DIR}/CircCoordinates"
CIRCTOOLS_CANONICAL_DCC_LINEAR_COUNT = f"{CIRCTOOLS_CANONICAL_DCC_DIR}/LinearCount"


# ---------------------------------------------------------------------------
# circtools helper functions (called from circtools.smk param lambdas).
# ---------------------------------------------------------------------------
def circtools_circtest_case_samples(comparison):
    """Case samples for a circtools circtest comparison, in declared order."""
    return list(CIRCTOOLS_CIRCTEST_COMPARISONS[comparison]["case"])


def circtools_circtest_control_samples(comparison):
    """Control samples for a circtools circtest comparison, in declared order."""
    return list(CIRCTOOLS_CIRCTEST_COMPARISONS[comparison]["control"])


def circtools_circtest_condition_names(comparison):
    """(case_group, control_group) tuple for -l on circtest/quickcheck."""
    comp = CIRCTOOLS_CIRCTEST_COMPARISONS[comparison]
    return (comp["case_group"], comp["control_group"])


def circtools_circtest_condition_columns(comparison):
    """1-based DCC sample column indices for circtest's -c flag.

    The per-comparison DCC subset has 3 metadata columns (Chr, Start, End)
    followed by case-then-control sample columns (see circtools.smk:
    circtools_dcc_subset rule). So the 1-based sample columns are simply
    4, 5, ..., 3+N where N is the total case+control count. We list ALL of
    them; the grouping is what tells circtest which is which.
    """
    n = (
        len(circtools_circtest_case_samples(comparison))
        + len(circtools_circtest_control_samples(comparison))
    )
    return ",".join(str(i) for i in range(4, 4 + n))


def circtools_circtest_grouping(comparison):
    """Comma-separated group ids for circtest's -g flag.

    Group 1 = case (matches the first name in condition_names), group 2 =
    control. Order MUST match circtools_circtest_condition_columns above and
    the sample order baked into the per-comparison DCC subset.
    """
    case_n = len(circtools_circtest_case_samples(comparison))
    control_n = len(circtools_circtest_control_samples(comparison))
    return ",".join(["1"] * case_n + ["2"] * control_n)


def circtools_circtest_num_replicates(comparison):
    """`-r` for circtest. The wrapper expects a single replicate count, but
    accepts imbalanced designs; using max(case, control) is the standard
    fix-up the upstream README also recommends.
    """
    return max(
        len(circtools_circtest_case_samples(comparison)),
        len(circtools_circtest_control_samples(comparison)),
    )


def circtools_source_star_log(sample):
    """Source Log.final.out path used by the quickcheck staging symlink rule.

    When the circtools detect arm is on, we have a dedicated STAR alignment
    that writes {sample}_Log.final.out under {OUTDIR}/circtools/star/{sample}/paired/.
    Otherwise, we reuse the CIRI3 STAR tree under {OUTDIR}/star/, which
    writes {sample}.Log.final.out (dot separator -- not directly compatible
    with quickcheck's lookup, which is why we stage a symlink renamed to a
    plain Log.final.out).
    """
    if CIRCTOOLS_DETECT_ENABLED:
        return f"{OUTDIR}/circtools/star/{sample}/paired/{sample}_Log.final.out"
    return f"{OUTDIR}/star/{sample}/{sample}.Log.final.out"


# ---------------------------------------------------------------------------
# Pre-computed target lists for the Snakefile's rule all.
# ---------------------------------------------------------------------------
CIRCTOOLS_DETECT_TARGETS = (
    [
        CIRCTOOLS_CANONICAL_DCC_CIRC_COUNT,
        CIRCTOOLS_CANONICAL_DCC_CIRC_COORD,
        CIRCTOOLS_CANONICAL_DCC_LINEAR_COUNT,
    ]
    if CIRCTOOLS_DETECT_ENABLED
    else []
)
CIRCTOOLS_CIRCTEST_TARGETS = (
    [
        f"{OUTDIR}/circtools/circtest/{comparison}/circtest.{ext}"
        for comparison in CIRCTOOLS_CIRCTEST_COMPARISON_NAMES
        for ext in ("csv", "xlsx", "pdf")
    ]
    if CIRCTOOLS_CIRCTEST_ENABLED
    else []
)
CIRCTOOLS_PRIMEX_TARGETS = (
    [f"{OUTDIR}/circtools/primex/primex.done"]
    if CIRCTOOLS_PRIMEX_ENABLED
    else []
)
CIRCTOOLS_QUICKCHECK_TARGETS = (
    [
        f"{OUTDIR}/circtools/quickcheck/{comparison}/quickcheck.pdf"
        for comparison in CIRCTOOLS_CIRCTEST_COMPARISON_NAMES
    ]
    if CIRCTOOLS_QUICKCHECK_ENABLED
    else []
)
CIRCTOOLS_INSTALL_TARGETS = (
    [CIRCTOOLS_INSTALL_MARKER]
    if CIRCTOOLS_ENABLED
    else []
)
