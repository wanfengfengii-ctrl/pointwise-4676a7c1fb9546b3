import argparse
import json
import sqlite3
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from pathlib import Path
from . import Archive, Cell, Target, RemapError


def json_value(value):
    if isinstance(value, Fraction):
        return str(value)
    if is_dataclass(value):
        return json_value(asdict(value))
    if isinstance(value, dict):
        return {k: json_value(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_value(v) for v in value]
    return value


def parser():
    p = argparse.ArgumentParser(description="Exact tetrahedral observation remapping")
    p.add_argument("--db", required=True)
    sub = p.add_subparsers(dest="command", required=True)
    s = sub.add_parser("source")
    s.add_argument("file")
    s = sub.add_parser("plan")
    s.add_argument("file")
    s = sub.add_parser("execute")
    s.add_argument("plan_id", type=int)
    s.add_argument("request_id")
    s = sub.add_parser("result")
    s.add_argument("result_id", type=int)
    s = sub.add_parser("head")
    s.add_argument("name")
    sub.add_parser("catalog")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        with Archive(args.db) as archive:
            if args.command == "source":
                data = json.loads(Path(args.file).read_text(encoding="utf-8"))
                cells = [Cell(c["id"], c["vertices"], c["density"])
                         for c in data["cells"]]
                result = {"source_revision": archive.put_source(data["name"], cells)}
            elif args.command == "plan":
                data = json.loads(Path(args.file).read_text(encoding="utf-8"))
                targets = [Target(t["id"], t["lower"], t["upper"])
                           for t in data["targets"]]
                result = {"plan_id": archive.save_plan(
                    data["name"], data["source_revision"], targets)}
            elif args.command == "execute":
                result = archive.execute_plan(args.plan_id, args.request_id)
            elif args.command == "result":
                result = archive.read_result(args.result_id)
            elif args.command == "head":
                result = {"source_revision": archive.source_head(args.name)}
            else:
                result = {"sources": archive.list_sources(), "plans": archive.list_plans(),
                          "results": archive.list_results()}
        print(json.dumps(json_value(result), ensure_ascii=False, sort_keys=True))
        return 0
    except (RemapError, sqlite3.Error, OSError, ValueError, KeyError, TypeError) as exc:
        print(json.dumps({"error": str(exc)}, ensure_ascii=False))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
