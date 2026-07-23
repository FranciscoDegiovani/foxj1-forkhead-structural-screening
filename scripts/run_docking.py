# Batch docking script using the standalone vina.exe binary (avoids the
# Boost compilation issue with the Python "vina" pip package on Windows).
# Runs multiple ligands in parallel across CPU cores, and checkpoints
# each result to CSV immediately so the run can be safely interrupted
# and resumed without losing progress.

import os
import csv
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

VINA_EXE = "vina.exe"
RECEPTOR = "data/processed/receptor/FOXJ1_receptor.pdbqt"
LIGAND_DIR = "data/processed/ligands_pdbqt"
OUT_DIR = "data/processed/docking_out"
RESULTS_CSV = "data/processed/docking_results.csv"

CENTER = ("-13.98", "3.32", "-4.04")
SIZE = ("22", "22", "22")
EXHAUSTIVENESS = "8"
CPU_PER_JOB = "2"          # cores dedicated to EACH parallel vina.exe run
N_PARALLEL_JOBS = 8        # how many ligands to dock simultaneously (16 cores / 2 per job)

os.makedirs(OUT_DIR, exist_ok=True)

def already_done():
    done = set()
    if os.path.exists(RESULTS_CSV):
        with open(RESULTS_CSV, newline="") as f:
            for row in csv.DictReader(f):
                done.add(row["chembl_id"])
    return done

def dock_one(chembl_id):
    ligand_path = f"{LIGAND_DIR}/{chembl_id}.pdbqt"
    out_path = f"{OUT_DIR}/{chembl_id}_out.pdbqt"
    log_path = f"{OUT_DIR}/{chembl_id}_log.txt"

    cmd = [
        VINA_EXE,
        "--receptor", RECEPTOR,
        "--ligand", ligand_path,
        "--center_x", CENTER[0], "--center_y", CENTER[1], "--center_z", CENTER[2],
        "--size_x", SIZE[0], "--size_y", SIZE[1], "--size_z", SIZE[2],
        "--exhaustiveness", EXHAUSTIVENESS,
        "--cpu", CPU_PER_JOB,
        "--out", out_path,
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        with open(log_path, "w") as f:
            f.write(result.stdout)
            f.write("\n--- STDERR ---\n")
            f.write(result.stderr)

        if result.returncode != 0:
            return (chembl_id, None, f"erro (codigo {result.returncode}): {result.stderr.strip()[:200]}")

        best_score = None
        for line in result.stdout.splitlines():
            if line.strip().startswith("1 "):  # first pose row in Vina's result table
                best_score = float(line.split()[1])
                break

        return (chembl_id, best_score, "ok")
    except Exception as e:
        return (chembl_id, None, f"erro: {e}")

def main():
    done = already_done()
    ligand_ids = sorted(
        f.replace(".pdbqt", "") for f in os.listdir(LIGAND_DIR) if f.endswith(".pdbqt")
    )
    pending = [lid for lid in ligand_ids if lid not in done]

    print(f"Total de ligantes: {len(ligand_ids)}")
    print(f"Ja processados: {len(done)}")
    print(f"Pendentes: {len(pending)}")

    write_header = not os.path.exists(RESULTS_CSV)
    with open(RESULTS_CSV, "a", newline="") as csvfile:
        writer = csv.writer(csvfile)
        if write_header:
            writer.writerow(["chembl_id", "best_affinity_kcal_mol", "status"])

        with ProcessPoolExecutor(max_workers=N_PARALLEL_JOBS) as executor:
            futures = {executor.submit(dock_one, lid): lid for lid in pending}
            completed = 0
            t0 = time.time()
            for future in as_completed(futures):
                chembl_id, score, status = future.result()
                writer.writerow([chembl_id, score, status])
                csvfile.flush()
                completed += 1
                if completed % 10 == 0:
                    elapsed = time.time() - t0
                    rate = completed / elapsed
                    remaining = (len(pending) - completed) / rate if rate > 0 else 0
                    print(f"{completed}/{len(pending)} concluidos "
                          f"({rate*60:.1f}/min) - ETA: {remaining/60:.0f} min")

    print("Docking em lote finalizado.")

if __name__ == "__main__":
    main()