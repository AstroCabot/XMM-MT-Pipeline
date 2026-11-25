import sys

from _common import OUT, REPRO, stage_dir

sys.path.insert(0, str(REPRO))

from analysis_05_radial_profile import data, run


def main():
    data.IMAGE = OUT / "04" / "image.fits"
    run.SPECTRUM = OUT / "03" / "result.json"
    run.MORPHOLOGY = OUT / "04" / "result.json"
    run.PPS = OUT / "02"
    run.OUTPUT = stage_dir(5)
    run.main()


if __name__ == "__main__":
    main()
