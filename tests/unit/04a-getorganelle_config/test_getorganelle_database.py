"""
test_getorganelle_database.py

Unit tests for the job-level GetOrganelle database check: that an index
built by another BLAST version (the real "Not a valid version 4
database" failure that silently sank a whole batch) is detected and
rebuilt with this environment's own tools, that missing databases are
fetched, and that a failure is reported as one -- all with a fake
runner, no real BLAST.
"""

# Imports
import importlib
import json
import logging
from pathlib import Path
import subprocess
import sys
import pytest
from mitopipeline.utils.assembly.getorganelle_database import GetOrganelleDatabase, default_config_dir

run_getorganelle_database = importlib.import_module("mitopipeline.exec.04a-getorganelle_config.run_getorganelle_database")

class _FakeTools:
    """
    A scripted environment: which tools exist, whether blastdbcmd can
    read the index (and whether makeblastdb fixes that), and what each
    command does to the filesystem.
    """
    def __init__(self, db: GetOrganelleDatabase | None = None, index_readable: bool = True,
                 makeblastdb_ok: bool = True, missing_tools: tuple[str, ...] = (),
                 download_ok: bool = True) -> None:
        self.db = db
        self.index_readable = index_readable
        self.makeblastdb_ok = makeblastdb_ok
        self.missing_tools = missing_tools
        self.download_ok = download_ok
        self.commands: list[list[str]] = []

    def which(self, tool: str) -> str | None:
        return None if tool in self.missing_tools else f"/env/bin/{tool}"

    def run(self, argv):
        argv = [str(a) for a in argv]
        self.commands.append(argv)
        tool = argv[0]
        if tool == "blastn":
            return subprocess.CompletedProcess(argv, 0, stdout = "blastn: 2.16.0+\n Package: blast 2.16.0\n", stderr = "")
        if tool == "blastdbcmd":
            return subprocess.CompletedProcess(argv, 0 if self.index_readable else 2, stdout = "",
                                               stderr = "" if self.index_readable else "BLAST Database error: Error: Not a valid version 4 database.")
        if tool == "makeblastdb":
            if self.makeblastdb_ok:
                self.index_readable = True
                return subprocess.CompletedProcess(argv, 0, stdout = "Building a new DB", stderr = "")
            return subprocess.CompletedProcess(argv, 1, stdout = "", stderr = "makeblastdb failed")
        if tool == "bowtie2-build":
            Path(argv[2] + ".1.bt2").write_bytes(b"")
            return subprocess.CompletedProcess(argv, 0, stdout = "", stderr = "")
        if tool == "get_organelle_config.py":
            if self.download_ok and self.db is not None:
                _make_fastas(self.db)
                self.index_readable = True
                Path(str(self.db.seed_prefix) + ".1.bt2").write_bytes(b"")
                return subprocess.CompletedProcess(argv, 0, stdout = "", stderr = "")
            return subprocess.CompletedProcess(argv, 1, stdout = "", stderr = "download failed")
        raise AssertionError(f"unexpected command {argv}")

    def ran(self, tool: str) -> list[list[str]]:
        return [c for c in self.commands if c[0] == tool]

def _make_fastas(db: GetOrganelleDatabase) -> None:
    db.label_fasta.parent.mkdir(parents = True, exist_ok = True)
    db.seed_fasta.parent.mkdir(parents = True, exist_ok = True)
    db.label_fasta.write_text(">label\nACGT\n")
    db.seed_fasta.write_text(">seed\nACGT\n")

def _database(tmp_path: Path, **tool_kwargs) -> tuple[GetOrganelleDatabase, _FakeTools]:
    tools = _FakeTools(**tool_kwargs)
    db = GetOrganelleDatabase(organelle_type = "animal_mt", logger = logging.getLogger("test"),
                              config_dir = tmp_path / ".GetOrganelle", runner = tools.run, which = tools.which)
    tools.db = db
    return db, tools

def _seed_index(db: GetOrganelleDatabase) -> None:
    db.seed_prefix.parent.mkdir(parents = True, exist_ok = True)
    Path(str(db.seed_prefix) + ".1.bt2").write_bytes(b"")

# ---- GetOrganelleDatabase.ensure() ----

def test_usable_databases_are_left_alone(tmp_path: Path) -> None:
    db, tools = _database(tmp_path)
    _make_fastas(db); _seed_index(db)

    result = db.ensure()
    assert result.success is True
    assert result.label_index_rebuilt is False
    assert result.seed_index_rebuilt is False
    assert result.downloaded is False
    assert result.blast_version == "blastn: 2.16.0+"
    assert tools.ran("makeblastdb") == []
    assert tools.ran("get_organelle_config.py") == []

def test_an_index_unreadable_by_this_blast_is_rebuilt_with_this_blast(tmp_path: Path) -> None:
    """Regression for the real failure: the LabelDatabase index under the
    shared ~/.GetOrganelle was built by another BLAST version, so
    blastdbcmd (and every slim_graph.py call inside GetOrganelle) died
    with "Not a valid version 4 database" and every sample ended in "No
    valid assembly graph found!". The check must notice and rebuild the
    index with the environment's own makeblastdb -- never leave it to a
    user to diagnose from a swallowed GetOrganelle error."""
    db, tools = _database(tmp_path, index_readable = False)
    _make_fastas(db); _seed_index(db)

    result = db.ensure()
    assert result.success is True
    assert result.label_index_rebuilt is True
    rebuilds = tools.ran("makeblastdb")
    assert len(rebuilds) == 1
    assert rebuilds[0][1:] == ["-in", str(db.label_fasta), "-dbtype", "nucl", "-out", str(db.label_prefix)]
    # Verified again after the rebuild, not assumed.
    assert len(tools.ran("blastdbcmd")) == 2
    assert "rebuilt: LabelDatabase BLAST index" in result.reason

def test_a_rebuild_that_still_is_not_readable_is_a_failure(tmp_path: Path) -> None:
    db, tools = _database(tmp_path, index_readable = False, makeblastdb_ok = False)
    _make_fastas(db); _seed_index(db)

    result = db.ensure()
    assert result.success is False
    assert "LabelDatabase" in result.reason

def test_missing_databases_are_downloaded_first(tmp_path: Path) -> None:
    db, tools = _database(tmp_path)

    result = db.ensure()
    assert result.success is True
    assert result.downloaded is True
    assert tools.ran("get_organelle_config.py") == [["get_organelle_config.py", "--add", "animal_mt"]]
    assert db.label_fasta.is_file() and db.seed_fasta.is_file()

def test_a_failed_download_is_a_failure(tmp_path: Path) -> None:
    db, _ = _database(tmp_path, download_ok = False)
    result = db.ensure()
    assert result.success is False
    assert "get_organelle_config.py" in result.reason

def test_a_missing_seed_index_is_built(tmp_path: Path) -> None:
    db, tools = _database(tmp_path)
    _make_fastas(db)  # no bowtie2 index

    result = db.ensure()
    assert result.success is True
    assert result.seed_index_rebuilt is True
    assert tools.ran("bowtie2-build") == [["bowtie2-build", str(db.seed_fasta), str(db.seed_prefix)]]

def test_missing_tools_are_reported_without_running_anything(tmp_path: Path) -> None:
    db, tools = _database(tmp_path, missing_tools = ("makeblastdb", "blastdbcmd"))
    _make_fastas(db); _seed_index(db)

    result = db.ensure()
    assert result.success is False
    assert "makeblastdb" in result.reason and "blastdbcmd" in result.reason
    assert tools.ran("get_organelle_config.py") == [] and tools.ran("makeblastdb") == []

def test_default_config_dir_honours_getorg_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.delenv("GETORG_PATH", raising = False)
    assert default_config_dir() == Path.home() / ".GetOrganelle"
    monkeypatch.setenv("GETORG_PATH", str(tmp_path / "custom"))
    assert default_config_dir() == tmp_path / "custom"

# ---- run_getorganelle_database.main() ----

def _argv(tmp_path: Path) -> list[str]:
    return ["run_getorganelle_database.py", "--organelle-type", "animal_mt",
            "--config-dir", str(tmp_path / ".GetOrganelle"),
            "--output-dir", str(tmp_path / "04a"),
            "--log-file", str(tmp_path / "logs" / "db.log"),
            "--global-log-file", str(tmp_path / "logs" / "global.log")]

def test_main_records_the_outcome_and_exits_zero_when_usable(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    class _FakeDatabase:
        def __init__(self, organelle_type, logger, config_dir): self.organelle_type, self.config_dir = organelle_type, config_dir
        def ensure(self):
            from mitopipeline.utils.assembly.getorganelle_database import GetOrganelleDatabaseResult
            return GetOrganelleDatabaseResult(success = True, config_dir = self.config_dir, organelle_type = self.organelle_type,
                                              blast_version = "blastn: 2.16.0+", label_index_rebuilt = True,
                                              reason = "databases ready; rebuilt: LabelDatabase BLAST index",
                                              actions = ["makeblastdb -in x -dbtype nucl -out y"])
    monkeypatch.setattr(run_getorganelle_database, "GetOrganelleDatabase", _FakeDatabase)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path))

    assert run_getorganelle_database.main() == 0
    record = json.loads((tmp_path / "04a" / "getorganelle_database.json").read_text())
    assert record["success"] is True
    assert record["label_index_rebuilt"] is True
    assert record["organelle_type"] == "animal_mt"
    assert record["actions"] == ["makeblastdb -in x -dbtype nucl -out y"]

def test_main_exits_nonzero_when_databases_are_unusable(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Mandatory, unlike MitoFish: an unusable LabelDatabase means no
    assembly can succeed, so the run stops here with the reason."""
    class _FakeDatabase:
        def __init__(self, organelle_type, logger, config_dir): self.organelle_type, self.config_dir = organelle_type, config_dir
        def ensure(self):
            from mitopipeline.utils.assembly.getorganelle_database import GetOrganelleDatabaseResult
            return GetOrganelleDatabaseResult(success = False, config_dir = self.config_dir, organelle_type = self.organelle_type,
                                              reason = "could not build a LabelDatabase BLAST index readable by this environment")
    monkeypatch.setattr(run_getorganelle_database, "GetOrganelleDatabase", _FakeDatabase)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path))

    assert run_getorganelle_database.main() == 1
    record = json.loads((tmp_path / "04a" / "getorganelle_database.json").read_text())
    assert record["success"] is False
    assert "LabelDatabase" in record["reason"]
