# workflow/rules/helpers.smk
import os

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


def maybe_temp(path):
    return path if KEEP_BAM else temp(path)


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

DEG_ACTIVE = DEG_RUN_DESEQ2 or CIRI3_DE_ENABLED

if DEG_ACTIVE:
    if len(DEG_GROUP_NAMES) < 2:
        raise ValueError("DEG/CIRI3 differential analysis is enabled, but fewer than 2 deg.groups are configured.")

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
                "biological replicates per group are required for enabled DEG/CIRI3 comparisons."
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

motif_cfg = config.get("motif", {})
MOTIF_ENABLED = _as_bool(motif_cfg.get("enabled", True), default=True)

if DEG_ACTIVE and not DEG_COMPARISON_NAMES:
    raise ValueError("DEG/CIRI3 differential analysis is enabled, but no comparisons could be generated.")

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
