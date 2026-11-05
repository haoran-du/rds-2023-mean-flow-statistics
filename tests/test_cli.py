import importlib
from pathlib import Path
import pytest
from rds2023.cli import main


@pytest.mark.parametrize("failed", [[], [{"name": "known legacy mismatch", "pass": False}]])
def test_validation_exit_status(tmp_path, monkeypatch, failed):
    source = tmp_path/"DATA"
    source.mkdir()
    module = importlib.import_module("rds2023.raw")
    monkeypatch.setattr(module, "validate_xz", lambda *a, **k: {"passed": 2, "failed": failed})
    args = ["validate-xz", "--data-root", str(source), "--output", str(tmp_path/"out")]
    if failed:
        with pytest.raises(SystemExit) as exc:
            main(args)
        assert exc.value.code == 1
    else:
        main(args)


def test_source_is_read_only(tmp_path):
    with pytest.raises(SystemExit) as exc:
        main(["validate-xz", "--data-root", str(tmp_path), "--output", str(tmp_path/"out")])
    assert exc.value.code == 2
    assert not (tmp_path/"out").exists()


def test_no_cross_paper_imports():
    package = importlib.import_module("rds2023")
    for path in Path(package.__file__).parent.glob("*.py"):
        source = path.read_text(encoding="utf-8")
        for forbidden in ("from rds.", "import rds", 'rds2024', 'rds2025'):
            assert forbidden not in source
