from pathlib import Path

from register_primers import discover_plasmid_files


def test_excludes_files_under_backup_dir(tmp_path):
    plasmid_dir = tmp_path / "PlasmidList"
    plasmid_dir.mkdir()
    real = plasmid_dir / "V001-01.example.dna"
    real.write_bytes(b"")

    backup_dir = plasmid_dir.parent / "_primer_registration_backups" / "20260101_000000"
    backup_dir.mkdir(parents=True)
    backup_copy = backup_dir / "V001-01.example.dna"
    backup_copy.write_bytes(b"")

    all_files = sorted(plasmid_dir.parent.rglob("*.dna"))
    assert backup_copy in all_files  # sanity check: rglob does see it

    discovered = discover_plasmid_files(all_files, backup_dir)
    assert discovered == [real]


def test_excludes_any_directory_named_backups_even_if_not_the_configured_backup_dir(tmp_path):
    """Regression test: an earlier run's leftover backup folder, sitting
    *inside* plasmid-dir from before this fix existed, must never be
    re-ingested as plasmid input on a later run -- even though it is not the
    *current* run's --backup-dir."""
    plasmid_dir = tmp_path / "PlasmidList"
    plasmid_dir.mkdir()
    real = plasmid_dir / "V001-01.example.dna"
    real.write_bytes(b"")

    stale_backup_dir = plasmid_dir / "_primer_registration_backups" / "20260101_000000"
    stale_backup_dir.mkdir(parents=True)
    (stale_backup_dir / "V001-01.example.dna").write_bytes(b"")

    all_files = sorted(plasmid_dir.rglob("*.dna"))
    assert len(all_files) == 2  # sanity check: both are on disk

    current_backup_dir = plasmid_dir.parent / "_primer_registration_backups" / "20260102_000000"
    discovered = discover_plasmid_files(all_files, current_backup_dir)
    assert discovered == [real]


def test_no_backup_dir_present_returns_everything(tmp_path):
    plasmid_dir = tmp_path / "PlasmidList"
    plasmid_dir.mkdir()
    real = plasmid_dir / "V001-01.example.dna"
    real.write_bytes(b"")

    all_files = sorted(plasmid_dir.rglob("*.dna"))
    backup_dir = plasmid_dir.parent / "_primer_registration_backups" / "20260101_000000"
    assert discover_plasmid_files(all_files, backup_dir) == [real]
