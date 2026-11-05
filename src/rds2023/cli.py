import argparse
import json
import sys
from pathlib import Path


def main(argv=None):
    parser = argparse.ArgumentParser(description="Paper-specific RdS reproduction; source data are read-only")
    sub = parser.add_subparsers(dest="command", required=True)
    repro = sub.add_parser("reproduce", help="Recover numerical figures 2-10")
    repro.add_argument("--normalization", choices=["paper", "legacy"], default="paper")
    raw = sub.add_parser("validate-xz", help="Recompute all XZ moments and fluxes")
    raw.add_argument("--source", choices=["mat", "txt"], default="mat")
    raw.add_argument("--batch", type=int, default=100)
    raw.add_argument("--cases", nargs="+", choices=["1b", "2b", "3b", "6b", "ca"])
    hwa = sub.add_parser("validate-hwa", help="Recompute HWA moments and spectra")
    hwa.add_argument("--cases", nargs="+", choices=["1b", "2b", "3b", "6b", "ca"])
    hwa.add_argument("--b6-policy", choices=["saved-zero", "recorded-replacement"], default="saved-zero")
    commands = (repro, raw, hwa)
    for command in commands:
        command.add_argument("--data-root", type=Path, required=True)
        command.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    root, output = args.data_root.resolve(strict=True), args.output.resolve()
    if output == root or root in output.parents:
        parser.error("Output must be outside the read-only source data tree")
    if getattr(args, "batch", 1) < 1:
        parser.error("--batch must be positive")
    if args.command == "reproduce":
        from .analysis import reproduce
        result = reproduce(args.data_root, args.output, args.normalization)
    elif args.command == "validate-xz":
        from .raw import validate_xz
        result = validate_xz(args.data_root, args.output, args.cases, args.batch, args.source)
    else:
        from .hwa import validate_hwa
        result = validate_hwa(args.data_root, args.output, args.cases, b6_policy=args.b6_policy)
    print(json.dumps(result, indent=2))
    if result.get("failed"):
        sys.exit(1)


if __name__ == "__main__":
    main()
