rule geo_download:
    output: "data/raw/{acc}_RAW.tar"
    params:
        url=lambda w: f"https://ftp.ncbi.nlm.nih.gov/geo/series/{w.acc[:-3]}nnn/{w.acc}/suppl/{w.acc}_RAW.tar"
    shell:
        """
        mkdir -p $(dirname {output})
        wget -q --show-progress -c -O {output} "{params.url}"
        """