{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# FNS Dogleg Duct Streaming Experiment Benchmark"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "In this folder the FNS Dogleg Duct Streaming Experiment is modeled in OpenMC for benchmark purposes. The experiment is part of the SINBAD database and, alongside experimental results, previous MCNP benchmark results are provided. The main objective is to benchmark OpenMC main capabilities in fusion relevant environments/models. More specifically, the present case study is a good candidate for testing and benchmarking OpenMC performances and capabilities with the implementation of the weight windows variance reduction technique."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "List of files:\n",
    "- ... (Python API OpenMC input providing geometry, materials, tallies and settings without implementation of variance reduction techniques)\n",
    "- ... (OpenMC input with weight window variance reduction technique implemented)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "The experiment consists of a Deuterium-Tritium (DT) fusion neutron source and a set of 11 detectors measuring reaction rates and neutron energy spectra. The neutron source is made by an deuteron acceletator and a Ti-T target assembly. The DT neutron source is placed at the center of a 296x296x450 cm room. An opening is present on one wall of the room. In the opening a 140x170x180 cm iron assembly is positioned. A doubly bent dogleg duct 30x30 cm in cross section was shaped through the assembly. The first duct horizontal leg of 115 cm is set at the same height of the neutron source. The second leg of 60 cm is connected vertically to the first with a right angle, and the third leg, long 65 cm, is horizontally connected to the second one. The detectors are placed in the duct and behind the streaming assembly. The source and the first detector are 170 cm apart. Detectors consist of a spectrometers and activation foils."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Bulding the model"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "The OpenMC model has been built mainly translating the MCNP input files provided in the SINBAD database. Another possible soultion would have been the application of the [openmc_mcnp-adapter](https://github.com/openmc-dev/openmc_mcnp_adapter) package. However, it does not provide the Python API input. The translation to a OpenMC Python API is quite straighforward for materials and geometry while the source and the tallies required a different approache."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Source"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "MCNP input files rely on a routine that models the entire accelerated deuterons - Ti-T target interaction specifically developed for MCNP. Other similar routines are available for MCNP (e.g. the Frascati Neutron Generator based one). However, a different modeling choice is necessary for OpenMC, as charged particles transport is not available. Currently, a\n",
    "<openmc.stats.Muir> point isotropic source is applied:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# source definition\n",
    "source = openmc.Source()\n",
    "source.particle = 'neutron'\n",
    "source.space = openmc.stats.Point([0,0,0])\n",
    "source.angle = openmc.stats.Isotropic()\n",
    "source.energy = openmc.stats.Muir()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Currently, such source performs relatively well."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "In the MCNP input the whole room and streaming assembly have been rotated around two axis in order to meet the source routine requirements. The current OpenMC input files perform the same rotations to match the MCNP input as much as possible. However, if the neutron source remains isotropic the rotations are not necessary and can be removed."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Tallies"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "In the experiment 4 cm diameter detectors filled with Ne-213 were used as spectrometer. Detectors were also equipped with Ni93, In155 and Au197 activation foils. Activation foils were applied to measure Nb93(n,2n)Nb92m, In115(n,n')In115m and Au197(n,g)Au198 reaction rates. Nevertheless, sizes and quantities of activation foils are not provided."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "MCNP input files model the detectors as 4x4x1 cm cells filled with air. Reaction rates rely on the Tally Multiplier \"fm\" feature of MCNP. Still, this choice have already been questioned in the \"THE QUALITY ASSESSMENT OF THE FNS BENCHMARK EXPERIMENTS\" document provided in the SINBAD database, along with the lack of information for consistently modeling the activation foils."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Currently, to get along with MCNP input files modeling choices we are filling the detector cells with the density modified activation foil materials. More specifically, the air density is applied:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# m33 Nb-93(n,2n)Nb-92m\n",
    "nb93 = openmc.Material(material_id=33, name='nb93')\n",
    "nb93.add_nuclide('Nb93', 1.0, 'ao')\n",
    "# nb93.set_density('g/cm3', 18.57)\n",
    "nb93.set_density('g/cm3', 1.2e-3)  # air density"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Then, reaction rate results are multiplied by a coefficient that scales the number of atoms present in the cell with the air density with the number of atoms present in a tipycal activation foil (that has a tipically a volume of the order of hundreds of cubic millimiters), or integers of that, assuming more than one activation foil could have been used. This brings to results close to the experimental ones usually. This is a guessing game already and it requires several runs of the same model (one per activation foil material)."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Alternatively it is possible to model a mixed material for all the activation foils, still with air density:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Mixed activation foil detector material\n",
    "mixed_detector = openmc.Material.mix_materials([nb93, in115, au197, air], [1/3, 1/3, 1/3, 0.00], 'vo')"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "And then calibrate all the reaction rate results on the first detector in the postprocessing phase."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Variance reduction"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "MCNP input files rely on cell importance sampling. Weight windows are instead applied in OpenMC."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "coming soon..."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3.8.10 64-bit",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.10"
  },
  "orig_nbformat": 4,
  "vscode": {
   "interpreter": {
    "hash": "916dbcbb3f70747c44a77c7bcd40155683ae19c65e1c03b4aa3499c5328201f1"
   }
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
