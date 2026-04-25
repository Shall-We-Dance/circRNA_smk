import argparse
import shutil
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Normalize CIRI3 per-sample outputs to workflow-declared paths."
    )
    parser.add_argument("--sample", required=True)
    parser.add_argument("--result", required=True)
    parser.add_argument("--bsj", required=True)
    parser.add_argument("--fsj", required=True)
    parser.add_argument("--tmpdir", required=True)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def log_message(log_path, message):
    with log_path.open("a") as fh:
        fh.write(f"[normalize_ciri3_outputs] {message}\n")


def existing_nonempty(paths):
    out = []
    seen = set()
    for path in paths:
        path = Path(path)
        if path in seen:
            continue
        seen.add(path)
        if path.is_file() and path.stat().st_size > 0:
            out.append(path)
    return out


def limited_tree(paths):
    lines = []
    seen = set()
    for root in paths:
        root = Path(root)
        if not root.exists() or root in seen:
            continue
        seen.add(root)
        if root.is_file():
            lines.append(f"{root}\t{root.stat().st_size}")
            continue
        for path in sorted(p for p in root.rglob("*") if p.is_file())[:200]:
            lines.append(f"{path}\t{path.stat().st_size}")
    return lines


def ciri3_completed(log_path):
    try:
        text = log_path.read_text(errors="replace")
    except OSError:
        return False
    return "Collation of circRNA completed" in text or "Program run time:" in text


def raise_missing_output(kind, log_path, search_roots):
    tree = "\n".join(limited_tree(search_roots)) or "(no files found)"
    message = (
        f"CIRI3 exited without a detectable {kind} output and the log does not "
        "contain a normal completion marker. Check the CIRI3 log above for the "
        f"original error. Files searched:\n{tree}"
    )
    log_message(log_path, message)
    raise RuntimeError(message)


def write_placeholder(path, header, log_path, reason):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\t".join(header) + "\n")
    log_message(log_path, f"Created placeholder {path}: {reason}")


def copy_candidate(candidate, dest, log_path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if candidate.resolve() != dest.resolve():
        shutil.copyfile(candidate, dest)
        log_message(log_path, f"Copied {candidate} -> {dest}")
    else:
        log_message(log_path, f"Using CIRI3 output {dest}")


def normalize_matrix_header(path, sample):
    lines = path.read_text().splitlines()
    if not lines:
        lines = ["circRNA_ID\t" + sample]
    else:
        header = lines[0].split("\t") if lines[0] else []
        if len(header) == 0:
            header = ["circRNA_ID", sample]
        elif len(header) == 1:
            header.append(sample)
        else:
            header[1] = sample
        lines[0] = "\t".join(header)
    path.write_text("\n".join(lines) + "\n")


def result_candidates(result, sample, search_roots):
    exact_names = {
        result.name,
        result.name + ".ciri3",
        sample + ".ciri3",
        sample + ".ciri3.ciri3",
    }
    candidates = [
        result,
        result.with_name(result.name + ".ciri3"),
        result.parent / f"{sample}.ciri3",
        result.parent / f"{sample}.ciri3.ciri3",
    ]
    for root in search_roots:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            name = path.name
            if name in exact_names:
                candidates.append(path)
            elif sample in name and name.endswith(".ciri3") and "Matrix" not in name:
                candidates.append(path)
    return existing_nonempty(candidates)


def matrix_candidates(dest, result, sample, kind, search_roots):
    candidates = [
        dest,
        result.with_name(result.name + f".{kind}_Matrix"),
        result.with_name(result.name + f".ciri3.{kind}_Matrix"),
        result.parent / f"{sample}.ciri3.{kind}_Matrix",
        result.parent / f"{sample}.{kind}_Matrix",
        result.parent / f"{kind}_Matrix.txt",
    ]
    generic_names = {f"{kind}_Matrix", f"{kind}_Matrix.txt"}
    for root in search_roots:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            name = path.name
            parent_names = {p.name for p in path.parents}
            sample_in_path = sample in name or any(sample in p for p in parent_names)
            if kind + "_Matrix" in name and sample_in_path:
                candidates.append(path)
            elif name in generic_names and sample_in_path:
                candidates.append(path)
    return existing_nonempty(candidates)


def main():
    args = parse_args()
    sample = args.sample
    result = Path(args.result)
    bsj = Path(args.bsj)
    fsj = Path(args.fsj)
    tmpdir = Path(args.tmpdir)
    log_path = Path(args.log)
    search_roots = [result.parent, tmpdir]
    completed = ciri3_completed(log_path)

    result_hits = result_candidates(result, sample, search_roots)
    if result_hits:
        copy_candidate(result_hits[0], result, log_path)
    else:
        if not completed:
            raise_missing_output("result table", log_path, search_roots)
        tree = "\n".join(limited_tree(search_roots))
        log_message(
            log_path,
            "No non-empty CIRI3 result table found. Files searched:\n" + tree,
        )
        write_placeholder(
            result,
            ["circRNA_ID", "gene_id", "strand"],
            log_path,
            "CIRI3 produced no result table; treating sample as zero detected circRNAs.",
        )

    for dest, kind in [(bsj, "BSJ"), (fsj, "FSJ")]:
        hits = matrix_candidates(dest, result, sample, kind, search_roots)
        if hits:
            copy_candidate(hits[0], dest, log_path)
        else:
            if not completed:
                raise_missing_output(f"{kind} matrix", log_path, search_roots)
            tree = "\n".join(limited_tree(search_roots))
            log_message(
                log_path,
                f"No non-empty {kind} matrix found. Files searched:\n" + tree,
            )
            write_placeholder(
                dest,
                ["circRNA_ID", sample],
                log_path,
                f"CIRI3 produced no {kind} matrix; treating sample as zero detected circRNAs.",
            )
        normalize_matrix_header(dest, sample)


if __name__ == "__main__":
    main()
