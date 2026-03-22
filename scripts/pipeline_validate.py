import gzip

from pipeline_config import BEXT_HEADER, GROWTH_STEPLOG_HEADER, SHEAR_G_DATA_HEADER


def open_text(path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def prefer_compressed(compressed_path, plain_path):
    if compressed_path.is_file():
        return compressed_path
    if plain_path.is_file():
        return plain_path
    return compressed_path


def parse_cell_packing(path):
    with open_text(path) as handle:
        header = handle.readline().split()
        if len(header) != 2:
            raise ValueError("missing or invalid packing header")
        count = int(header[0])
        phi = float(header[1])
        if count <= 0:
            raise ValueError("packing must contain at least one cell")
        rows = 0
        for line in handle:
            text = line.strip()
            if not text:
                continue
            fields = text.split()
            if len(fields) != 5:
                raise ValueError(f"expected 5 columns per cell row, found {len(fields)}")
            for field in fields:
                float(field)
            rows += 1
    if rows != count:
        raise ValueError(f"header count {count} does not match {rows} cell rows")
    return {"count": count, "phi": phi}


def parse_lobe_packing(path):
    with open_text(path) as handle:
        frames = 0
        while True:
            header_line = handle.readline()
            if header_line == "":
                break
            if not header_line.strip():
                continue
            header = header_line.split()
            if len(header) != 2:
                raise ValueError("missing or invalid lobe-packing header")
            count = int(header[0])
            float(header[1])
            if count <= 0:
                raise ValueError("lobe packing must contain at least one particle")
            for _ in range(count):
                row = handle.readline()
                if row == "":
                    raise ValueError("truncated lobe-packing frame")
                fields = row.split()
                if len(fields) != 4:
                    raise ValueError(f"expected 4 columns per lobe row, found {len(fields)}")
                float(fields[0])
                float(fields[1])
                float(fields[2])
                int(fields[3])
            frames += 1
    if frames == 0:
        raise ValueError("lobe packing contains no frames")
    return frames


def parse_lineage_snapshot(path):
    with open_text(path) as handle:
        header = handle.readline()
        if not header.startswith("# index cell_id parent_id"):
            raise ValueError("unexpected lineage snapshot header")
        rows = 0
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 8:
                raise ValueError(f"expected 8 columns per lineage row, found {len(fields)}")
            int(fields[0])
            int(fields[1])
            int(fields[2])
            float(fields[3])
            float(fields[4])
            int(fields[5])
            int(fields[6])
            int(fields[7])
            rows += 1
    if rows <= 0:
        raise ValueError("lineage snapshot contains no rows")
    return rows


def parse_growth_stats(path, expected_rows):
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().split()
        if len(header) != 3:
            raise ValueError("missing or invalid stats header")
        row_count = int(header[0])
        float(header[1])
        float(header[2])
        if row_count != expected_rows:
            raise ValueError(f"stats header count {row_count} does not match expected {expected_rows}")
        rows = 0
        for line in handle:
            text = line.strip()
            if not text:
                continue
            fields = text.split()
            if len(fields) != 8:
                raise ValueError(f"expected 8 columns per stats row, found {len(fields)}")
            float(fields[0])
            float(fields[1])
            float(fields[2])
            int(fields[3])
            float(fields[4])
            float(fields[5])
            float(fields[6])
            float(fields[7])
            rows += 1
    if rows != expected_rows:
        raise ValueError(f"stats file has {rows} rows, expected {expected_rows}")


def parse_nc(path, expected_jamm_count, expected_frame_count):
    fields = path.read_text(encoding="utf-8").split()
    if len(fields) != 26:
        raise ValueError(f"expected 26 NC values, found {len(fields)}")
    values = list(map(float, fields))
    for offset, expected_count in ((0, expected_jamm_count), (13, expected_frame_count)):
        count = int(values[offset])
        if count != expected_count:
            raise ValueError(f"NC count {count} does not match expected {expected_count}")
        for index in range(offset, offset + 8):
            int(values[index])
        for index in range(offset + 8, offset + 13):
            float(values[index])


def parse_steplog(path, expected_final_count):
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip()
        if header != GROWTH_STEPLOG_HEADER:
            raise ValueError("unexpected STEPLOG header")
        rows = 0
        last_count = None
        saw_postjam = False
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            fields = text.split()
            if len(fields) != 10:
                raise ValueError(f"expected 10 columns per steplog row, found {len(fields)}")
            int(fields[0])
            last_count = int(fields[1])
            float(fields[2])
            float(fields[3])
            float(fields[4])
            float(fields[5])
            float(fields[6])
            before_jamming = int(fields[7])
            at_jamming = int(fields[8])
            above_jamming = int(fields[9])
            if before_jamming not in (0, 1) or at_jamming not in (0, 1) or above_jamming not in (0, 1):
                raise ValueError("invalid jamming-state flags in steplog")
            saw_postjam = saw_postjam or above_jamming == 1
            rows += 1
    if rows == 0:
        raise ValueError("STEPLOG contains no data rows")
    if last_count != expected_final_count:
        raise ValueError(f"STEPLOG final N {last_count} does not match expected {expected_final_count}")
    if not saw_postjam:
        raise ValueError("STEPLOG never reaches the post-jamming regime")


def parse_divlog(path):
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        if not header.startswith("# step phi parent_cell_id"):
            raise ValueError("unexpected DIVLOG header")
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 6:
                raise ValueError(f"expected 6 columns per DIVLOG row, found {len(fields)}")
            int(fields[0])
            float(fields[1])
            int(fields[2])
            int(fields[3])
            int(fields[4])
            int(fields[5])


def parse_transitions(path):
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        if not header.startswith("# step phi cell_id parent_id"):
            raise ValueError("unexpected TRANSITIONS header")
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 12:
                raise ValueError(f"expected 12 columns per TRANSITIONS row, found {len(fields)}")
            int(fields[0])
            float(fields[1])
            int(fields[2])
            int(fields[3])
            float(fields[4])
            float(fields[5])
            float(fields[6])
            int(fields[7])
            int(fields[8])
            int(fields[9])
            int(fields[10])


def parse_postjamm_summary(path):
    rows = 0
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        if not header.startswith("# step phi P N Nc Nf Nu Ziso total_growthrate chi_c"):
            raise ValueError("unexpected POSTJAMM_SUMMARY header")
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 19:
                raise ValueError(f"expected 19 columns per POSTJAMM_SUMMARY row, found {len(fields)}")
            int(fields[0])
            float(fields[1])
            float(fields[2])
            int(fields[3])
            int(fields[4])
            int(fields[5])
            int(fields[6])
            int(fields[7])
            float(fields[8])
            float(fields[9])
            int(fields[10])
            int(fields[11])
            int(fields[12])
            int(fields[13])
            int(fields[14])
            int(fields[15])
            int(fields[16])
            int(fields[17])
            int(fields[18])
            rows += 1
    if rows <= 0:
        raise ValueError("POSTJAMM_SUMMARY contains no rows")


def parse_initial_free_cohort(path):
    with open_text(path) as handle:
        header = handle.readline()
        if not header.startswith("# jam_cell_id jam_parent_id initial_bud_diameter"):
            raise ValueError("unexpected COHORT_INITIAL_FREE header")
        seen_cell_ids = set()
        for line in handle:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 21:
                raise ValueError(f"expected 21 columns per cohort row, found {len(fields)}")
            jam_cell_id = int(fields[0])
            if jam_cell_id in seen_cell_ids:
                raise ValueError("duplicate jam_cell_id in cohort file")
            seen_cell_ids.add(jam_cell_id)
            int(fields[1])
            initial_bud_diameter = float(fields[2])
            if initial_bud_diameter <= 0.0:
                raise ValueError("invalid initial_bud_diameter in cohort file")
            status_code = int(fields[3])
            if status_code not in (0, 1, 2):
                raise ValueError("invalid cohort status code")
            had_first_contact = int(fields[4])
            if had_first_contact not in (0, 1):
                raise ValueError("invalid had_first_contact flag in cohort file")
            first_contact_step = int(fields[5])
            first_contact_phi = float(fields[6])
            first_contact_delta = float(fields[7])
            completion_step = int(fields[8])
            completion_phi = float(fields[9])
            completion_delta = float(fields[10])
            division_step = int(fields[11])
            division_phi = float(fields[12])
            division_delta = float(fields[13])
            final_step = int(fields[14])
            final_phi = float(fields[15])
            final_bud_diameter = float(fields[16])
            final_delta = float(fields[17])
            final_bud_contacts = int(fields[18])
            final_bud_unconstrained = int(fields[19])
            final_bud_compressed = int(fields[20])

            if first_contact_delta < 0.0 or completion_delta < 0.0 or division_delta < 0.0:
                raise ValueError("negative cohort delta_bud_area")
            if final_delta < 0.0:
                raise ValueError("negative final_delta_bud_area")
            if final_bud_contacts < 0:
                raise ValueError("negative final_bud_contacts")
            if final_bud_unconstrained not in (0, 1):
                raise ValueError("invalid final_bud_unconstrained flag")
            if final_bud_compressed not in (0, 1):
                raise ValueError("invalid final_bud_compressed flag")

            if had_first_contact == 0 and (
                first_contact_step != 0 or first_contact_phi != 0.0 or first_contact_delta != 0.0
            ):
                raise ValueError("missing first-contact flag for populated first-contact fields")
            if had_first_contact == 1 and first_contact_step <= 0:
                raise ValueError("invalid first-contact step in cohort file")

            if status_code == 0:
                if completion_step != 0 or completion_phi != 0.0 or completion_delta != 0.0:
                    raise ValueError("active cohort row has populated completion fields")
                if division_step != 0 or division_phi != 0.0 or division_delta != 0.0:
                    raise ValueError("active cohort row has populated division fields")
                if final_step <= 0 or final_phi <= 0.0:
                    raise ValueError("active cohort row is missing final censoring state")
                if final_bud_diameter <= 0.0:
                    raise ValueError("active cohort row has invalid final_bud_diameter")
                if final_bud_unconstrained != 1:
                    raise ValueError("active cohort row must remain bud_unconstrained at final dphi")
            elif status_code == 1:
                if completion_step <= 0:
                    raise ValueError("completed cohort row is missing completion_step")
                if division_step != 0 or division_phi != 0.0 or division_delta != 0.0:
                    raise ValueError("completed cohort row has populated division fields")
                if any(value != 0 for value in (final_step, final_bud_contacts, final_bud_unconstrained, final_bud_compressed)):
                    raise ValueError("completed cohort row has populated final-state integer fields")
                if final_phi != 0.0 or final_bud_diameter != 0.0 or final_delta != 0.0:
                    raise ValueError("completed cohort row has populated final-state float fields")
            else:
                if division_step <= 0:
                    raise ValueError("divided cohort row is missing division_step")
                if completion_step != 0 or completion_phi != 0.0 or completion_delta != 0.0:
                    raise ValueError("divided cohort row has populated completion fields")
                if any(value != 0 for value in (final_step, final_bud_contacts, final_bud_unconstrained, final_bud_compressed)):
                    raise ValueError("divided cohort row has populated final-state integer fields")
                if final_phi != 0.0 or final_bud_diameter != 0.0 or final_delta != 0.0:
                    raise ValueError("divided cohort row has populated final-state float fields")


def parse_stdout_log(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty log file: {path.name}")


def parse_shear_g_data(path):
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip()
        if header != SHEAR_G_DATA_HEADER:
            raise ValueError("unexpected G_data header")
        rows = 0
        first_n = None
        previous_strain = None
        for line in handle:
            text = line.strip()
            if not text:
                continue
            fields = text.split()
            if len(fields) != 12:
                raise ValueError(f"expected 12 columns per G_data row, found {len(fields)}")
            strain = float(fields[0])
            float(fields[1])
            float(fields[2])
            count = int(fields[3])
            int(fields[4])
            int(fields[5])
            int(fields[6])
            int(fields[7])
            float(fields[8])
            int(fields[9])
            int(fields[10])
            int(fields[11])
            if count <= 0:
                raise ValueError("invalid particle count in G_data")
            if previous_strain is not None and strain < previous_strain:
                raise ValueError("non-monotonic strain sequence in G_data")
            previous_strain = strain
            if first_n is None:
                first_n = count
            elif count != first_n:
                raise ValueError("particle count changes across G_data rows")
            rows += 1
    if rows <= 0:
        raise ValueError("G_data contains no rows")
    return first_n


def parse_bext_data(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty file: {path.name}")
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip()
        if header != BEXT_HEADER:
            raise ValueError("unexpected B_ext header")
        rows = 0
        for line in handle:
            text = line.strip()
            if not text:
                continue
            fields = text.split()
            if len(fields) != 20:
                raise ValueError(f"expected 20 columns, found {len(fields)}")
            for index in (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 18, 19):
                float(fields[index])
            for index in (10, 11, 12, 13, 14, 15, 16, 17):
                int(fields[index])
            rows += 1
    if rows != 1:
        raise ValueError(f"expected exactly 1 B_ext data row, found {rows}")


def growth_done(paths):
    try:
        jamm_info = parse_cell_packing(prefer_compressed(paths["jamm_gz"], paths["jamm"]))
        frame_info = parse_cell_packing(prefer_compressed(paths["frame_gz"], paths["frame"]))
        if parse_lineage_snapshot(prefer_compressed(paths["lineage_jamm_gz"], paths["lineage_jamm"])) != jamm_info["count"]:
            raise ValueError("jammed lineage row count mismatch")
        if parse_lineage_snapshot(prefer_compressed(paths["lineage_frame_gz"], paths["lineage_frame"])) != frame_info["count"]:
            raise ValueError("final lineage row count mismatch")
        parse_growth_stats(paths["stats_jamm"], jamm_info["count"] * 2)
        parse_growth_stats(paths["stats_frame"], frame_info["count"] * 2)
        parse_nc(paths["nc"], jamm_info["count"], frame_info["count"])
        parse_steplog(paths["steplog"], frame_info["count"])
        parse_divlog(paths["divlog"])
        parse_transitions(paths["transitions"])
        parse_postjamm_summary(paths["postjamm_summary"])
        parse_initial_free_cohort(prefer_compressed(paths["cohort_gz"], paths["cohort"]))
        parse_stdout_log(paths["log"])
        if paths["traj_gz"].is_file():
            parse_lobe_packing(paths["traj_gz"])
        elif paths["traj"].is_file():
            parse_lobe_packing(paths["traj"])
        return True
    except (FileNotFoundError, OSError, ValueError):
        return False


def shear_done(paths):
    try:
        expected_count = parse_shear_g_data(paths["g_data"])
        parse_stdout_log(paths["log"])
        if paths["traj_gz"].is_file():
            frame_count = parse_lobe_packing(paths["traj_gz"])
            if frame_count <= 0:
                raise ValueError("invalid shear trajectory")
        elif paths["traj"].is_file():
            frame_count = parse_lobe_packing(paths["traj"])
            if frame_count <= 0:
                raise ValueError("invalid shear trajectory")
        return expected_count > 0
    except (FileNotFoundError, OSError, ValueError):
        return False


def bext_done(paths):
    try:
        parse_bext_data(paths["data"])
        parse_stdout_log(paths["log"])
        if paths["frame_gz"].is_file():
            parse_cell_packing(paths["frame_gz"])
        elif paths["frame"].is_file():
            parse_cell_packing(paths["frame"])
        return True
    except (FileNotFoundError, OSError, ValueError):
        return False
