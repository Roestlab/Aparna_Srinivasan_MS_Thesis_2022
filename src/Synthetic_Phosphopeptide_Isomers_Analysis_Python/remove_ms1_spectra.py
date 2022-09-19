
from pyopenms import *
import click
@click.command()
@click.option('--filename')
def main(filename):
    """Remove spectra with mz target = 600 and mz target = 800 with 200 Da wide windows from the Bruker Compass Data Analysis converted MzML file. """
    exp = MSExperiment()

    MzMLFile().load(filename, exp)

    spec = []
    for s in exp.getSpectra():
        if s.getPrecursors()[0].getMZ() != 600.0 and s.getPrecursors()[0].getMZ() != 800.0:
            spec.append(s)

    exp.setSpectra(spec)

    MzMLFile().store(f"{filename}_filtered.mzML", exp) 


if __name__ == "__main__":
    main()