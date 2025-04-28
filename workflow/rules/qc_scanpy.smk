rule qc_h5ad:
    input:
        tar="data/raw/{acc}_RAW.tar"
    output:
        "data/processed/{acc}.h5ad"
    run:
        import tarfile, pathlib, shutil, glob, os
        wd = pathlib.Path(f"data/raw/{wildcards.acc}")
        if wd.exists(): shutil.rmtree(wd)
        wd.mkdir(parents=True)
        tarfile.open(input.tar).extractall(wd)

        mtx = glob.glob(str(wd / "**/*matrix.mtx*"), recursive=True)
        if mtx:  # 10x branch
            mtx = mtx[0]
            pref = os.path.basename(mtx).split("_matrix")[0]
            genes = (glob.glob(str(wd/f"**/{pref}_features.tsv*"),recursive=True)
                     or glob.glob(str(wd/f"**/{pref}_genes.tsv*"),recursive=True))[0]
            bcs = glob.glob(str(wd/f"**/{pref}_barcodes.tsv*"),recursive=True)[0]
            shell("python scripts/qc_to_h5ad.py {mtx} {genes} {bcs} {output}")
        else:  # Counts-csv branch
            csvs = glob.glob(str(wd/"**/*Counts.csv.gz"),recursive=True)
            if not csvs:
                raise ValueError(f"Neither 10x nor *_Counts csv in {wd}")
            shell("python scripts/csv_to_h5ad.py {csvs[0]} {output}")