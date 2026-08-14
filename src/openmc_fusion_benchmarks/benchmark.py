from __future__ import annotations

import yaml
import h5py
from pathlib import Path
import warnings
from abc import ABC, abstractmethod
import numpy as np
import openmc
from cad_to_dagmc import CadToDagmc

from .validate_spec import validate_benchmark
from .benchmark_results import BenchmarkResults
from .report import ReportConfig, ResultSource, build_report, render_pdf, render_yaml
from .backends.openmc.tallies import (
    make_default_openmc_normalizer,
    save_openmc_statepoint_tallies,
)
from .uq.tmc_engine import tmc_engine


BENCHMARK_DIR = Path(__file__).parent / "benchmarks"
LFS_DIR = Path(__file__).parents[2] / "lfs"


class Benchmark(ABC):
    def __init__(self, name: str):
        self.name = name

        # Validate the benchmark specification
        validate_benchmark(name)

        with (BENCHMARK_DIR / f"{self.name}/specifications.yaml").open("r") as f:
            self._benchmark_spec = yaml.safe_load(f)

        self._read_metadata()

    @abstractmethod
    def _build_materials(self):
        """Build materials for the benchmark."""
        pass

    @abstractmethod
    def _build_geometry(self):
        """Build geometry for the benchmark."""
        pass

    @abstractmethod
    def _build_source(self):
        """Build source for the benchmark."""
        pass

    @abstractmethod
    def _build_settings(self):
        """Build settings for the benchmark."""
        pass

    @abstractmethod
    def _build_tallies(self):
        """Build tallies for the benchmark."""
        pass

    @abstractmethod
    def _build_model(self):
        """Build the whole model for the benchmark."""
        pass

    @abstractmethod
    def _postprocess(self):
        """Post-process the model after running."""
        pass

    @abstractmethod
    def run(
        self,
        *args,
        generate_report: bool = False,
        report_config: ReportConfig | None = None,
        **kwargs,
    ):
        """Run the benchmark simulation."""
        pass

    @abstractmethod
    def _uncertainty_quantification(self):
        """Perform uncertainty quantification for the benchmark."""
        pass

    def _read_metadata(self):
        """Read metadata from the benchmark specification."""
        metadata = self._benchmark_spec['metadata']

        lines = []
        lines.append(f"📘 Title: {metadata.get('title', 'N/A')}")
        lines.append("")
        lines.append(f"🔖 Type: {metadata.get('type', 'N/A')}")
        lines.append("")
        lines.append(f"📂 Category: {metadata.get('category', 'N/A')}")
        lines.append("")
        lines.append(f"🧮 Version: {metadata.get('version', 'N/A')}")
        lines.append("")
        lines.append(f"📝 Description: {metadata.get('description', 'N/A')}")
        lines.append(f"📅 Date: {metadata.get('date', 'N/A')}")

        location = metadata.get("location", {})
        lines.append("📍 Location:")
        lines.append(f"   - Facility: {location.get('facility', 'N/A')}")
        lines.append(f"   - City: {location.get('city', 'N/A')}")
        lines.append(f"   - Country: {location.get('country', 'N/A')}")
        lines.append("")

        references = metadata.get("references", [])
        lines.append("🔗 References:")
        for ref in references:
            lines.append(f"   - Title: {ref.get('title', 'N/A')}")
            if 'doi' in ref:
                lines.append(f"     DOI: {ref['doi']}")
            if 'url' in ref:
                lines.append(f"     URL: {ref['url']}")

        authors = metadata.get("authors", [])
        if authors:
            lines.append("")
            lines.append("👥 Authors:")
            for author in authors:
                name = author.get("name", "N/A")
                affiliation = author.get("affiliation", "N/A")
                email = author.get("email", "N/A")
                lines.append(f"   - {name} ({affiliation}, {email})")

        self._metadata = "\n".join(lines)

    @property
    def metadata(self):
        """Show metadata for the benchmark."""

        return self._metadata

    def _write_spec_snapshot(self, filename: str = "benchmark_results.h5") -> None:
        """Persist a snapshot of specifications.yaml into the results file."""
        path = Path(filename)
        if not path.exists():
            warnings.warn(
                f"Results file '{filename}' not found. Skipping spec snapshot.",
                UserWarning,
            )
            return

        spec_yaml = yaml.safe_dump(self._benchmark_spec, sort_keys=False)
        spec_bytes = spec_yaml.encode("utf-8")

        with h5py.File(path, "a") as handle:
            if "specifications" in handle:
                del handle["specifications"]
            group = handle.create_group("specifications")
            group.attrs["format"] = "yaml"
            group.attrs["benchmark_name"] = self.name
            group.create_dataset("yaml", data=np.bytes_(spec_bytes))

    def _write_run_metadata(
        self,
        code_name: str,
        code_version: str,
        nuclear_data_name: str | None = None,
        nuclear_data_version: str | None = None,
        geometry: str | None = None,
        filename: str = "benchmark_results.h5",
    ) -> None:
        """Persist run metadata into the results file."""
        path = Path(filename)
        if not path.exists():
            warnings.warn(
                f"Results file '{filename}' not found. Skipping run metadata.",
                UserWarning,
            )
            return

        with h5py.File(path, "a") as handle:
            if "run_metadata" in handle:
                del handle["run_metadata"]
            group = handle.create_group("run_metadata")
            group.attrs["code_name"] = str(code_name)
            group.attrs["code_version"] = str(code_version)
            if nuclear_data_name is not None:
                group.attrs["nuclear_data_name"] = str(nuclear_data_name)
            if nuclear_data_version is not None:
                group.attrs["nuclear_data_version"] = str(nuclear_data_version)
            if geometry is not None:
                group.attrs["geometry"] = str(geometry)

    def _generate_report(self, report_config: ReportConfig | None = None) -> None:
        """Generate a report from the current benchmark results file."""
        results_path = Path("benchmark_results.h5")
        if not results_path.exists():
            warnings.warn("benchmark_results.h5 not found. Skipping report generation.", UserWarning)
            return

        print("Generating validation report...")

        if report_config is None:
            report_config = ReportConfig(output_dir=Path("report"), include_yaml=True, include_pdf=True)

        sources: list[ResultSource] = []
        calculation = BenchmarkResults.from_file(results_path)
        sources.append(ResultSource(name="calculation", kind="calculation", results=calculation))

        try:
            reference = BenchmarkResults.from_database(self.name, filename="experiment.h5")
            sources.append(ResultSource(name="reference", kind="experiment", results=reference))
        except Exception:
            pass

        report = build_report(sources, report_config)
        output_dir = Path(report_config.output_dir)
        plots_dir = output_dir / "plots"
        if report_config.include_yaml:
            render_yaml(report, output_dir / "report.yaml")
        if report_config.include_pdf:
            render_pdf(report, output_dir / "report.pdf", plots_dir)


class OpenmcBenchmark(Benchmark):
    def __init__(self, name: str):
        super().__init__(name)
        self._materials = None
        self._geometry = None
        self._settings = None
        self._tallies = None

        self.model = self._build_model()

    def _build_materials(self):
        # Implement the logic to build materials for OpenMC
        material_data = self._benchmark_spec['materials']

        fraction_map = {'atomic': 'ao', 'weight': 'wo'}

        materials = openmc.Materials()
        for m in material_data:
            mat = openmc.Material(material_id=m['id'], name=m['name'])
            mat.set_density(m['density']['units'], m['density']['value'])

            # Ensure fraction type is valid
            fraction_type = m['composition']['fraction_type']
            if fraction_type not in fraction_map:
                raise ValueError(f"Invalid fraction type: {fraction_type}")

            ft = fraction_map[fraction_type]
            add_method = mat.add_element if m['composition']['composition_type'] == 'element' else mat.add_nuclide

            for k, v in m['composition']['data'].items():
                add_method(k, v, ft)

            materials.append(mat)

        return materials

    def _build_geometry(self):
        def build_mesh(cad_file: str, material_tags, set_size: dict, global_mesh_size_min: float, global_mesh_size_max: float, mesh_file: str = "mesh.h5m"):

            # Instantiate the CadToDagmc model
            model = CadToDagmc()
            # Load the STEP file and assign material tags
            model.add_stp_file(filename=cad_file,
                               material_tags=material_tags, scale_factor=.1)

            # Generate the mesh
            model.export_dagmc_h5m_file(imprint=True, min_mesh_size=global_mesh_size_min,
                                        max_mesh_size=global_mesh_size_max, set_size=set_size, filename=mesh_file)

        # Implement the logic to build geometry for OpenMC
        geometry_data = self._benchmark_spec['geometry']

        # Preprocess for mesh generation
        meshing = geometry_data['meshing']
        # Sort volumes by rising id
        volumes = sorted(meshing['volumes'], key=lambda x: x['id'])
        # Get material tags out of sorted volumes
        material_tags = [entry['material'] for entry in volumes]
        # Build the set size: specific mesh size for some volumes
        set_size = {v['id']: v['mesh_size']
                    for v in volumes if 'mesh_size' in v}
        # Get global mesh sizes
        global_mesh_size_min = meshing['global_mesh_size_min']
        global_mesh_size_max = meshing['global_mesh_size_max']

        # Get the STEP file
        cad_file = LFS_DIR / "benchmarks" / f"{geometry_data['cad_file']}"

        # Generate the mesh if mesh.h5m not already present
        if Path("mesh.h5m").exists():
            warnings.warn(
                f"Mesh file already exists. Skipping mesh generation.")
        else:
            build_mesh(cad_file=cad_file, material_tags=material_tags, set_size=set_size,
                       global_mesh_size_min=global_mesh_size_min, global_mesh_size_max=global_mesh_size_max)

        # download the h5m file
        # download_from_drive(benchmark_name=self.name, file_format='h5m')
        mesh_file = Path("mesh.h5m")
        # Implement the logic to build geometry for OpenMC
        dag_universe = openmc.DAGMCUniverse(
            mesh_file).bounded_universe(starting_id=90000)

        return openmc.Geometry(root=dag_universe)

    def _build_source(self):
        
        source_data = self._benchmark_spec['sources']

        def energy_conversion(values, units):
            values = np.array(values)
            if units == 'eV':
                return values
            elif units == 'keV':
                return values * 1e3
            elif units == 'MeV':
                return values * 1e6
            elif units == 'GeV':
                return values * 1e9
            else:
                raise ValueError(f"Unsupported energy unit: {units}")

        def angular_conversion(values, units):
            values = np.array(values)
            if units == 'degrees':
                return values * np.pi / 180
            elif units == 'radians':
                return values
            else:
                raise ValueError(f"Unsupported angle unit: {units}")

        sources = []
        # More than one source is possible
        for source in source_data:
            # Handle source particle type
            particle = source['particle']

            # Handle source spatial distribution
            if source['spatial_distribution']['type'] == 'point':
                center = source['spatial_distribution']['center']
                space = openmc.stats.Point(center)
            elif source['spatial_distribution']['type'] == 'box':
                raise NotImplementedError(
                    'Box source distribution not implemented yet.')
            elif source['spatial_distribution']['type'] == 'sphere':
                raise NotImplementedError(
                    'Sphere source distribution not implemented yet.')
            elif source['spatial_distribution']['type'] == 'cylinder':
                raise NotImplementedError(
                    'Cylinder source distribution not implemented yet.')
            elif source['spatial_distribution']['type'] == 'cartesian':
                raise NotImplementedError(
                    'Cartesian source distribution not implemented yet.')
            elif source['spatial_distribution']['type'] == 'mesh':
                raise NotImplementedError(
                    'Mesh source distribution not implemented yet.')
            else:
                raise ValueError(
                    f"Unsupported spatial distribution type: {source['spatial_distribution']['type']}")

            # # Handle if source is associated with a domain (e.g. a cell, volume or material)
            # if source['spatial_distribution']['domain'] is not None:
            #     raise NotImplementedError(
            #         'Source domain not implemented yet.')

            # Handle angular and energy distributions
            angular_sources = []
            # Openmc needs one source per angle bin:
            angles = source['angular_energy_distribution']
            abins = np.cos(angles['angle']['bins'])
            for i in range(len(abins)-1):
                lb = angular_conversion(
                    abins[i], angles['angle']['units'])
                ub = angular_conversion(
                    abins[i+1], angles['angle']['units'])
                mu = openmc.stats.Uniform(a=lb, b=ub)  # polar angle
                # Azimuthal angle --> STILL TO HANDLE
                phi = openmc.stats.Uniform(0, 2*np.pi)
                reference_uvw = angles['polar_direction']
                # Polar-azimuthal distribution
                angle = openmc.stats.PolarAzimuthal(mu, phi, reference_uvw)
                # Energy distribution
                evalues = energy_conversion(
                    angles['energy']['values'], angles['energy']['units'])
                # Weights
                weights = np.array(angles['weights'])

                # Check weights has the right shape
                # Weights must be a 2D array (rows, columns)
                if weights.ndim != 2:
                    raise ValueError(
                        'Weights must be a 2D array (rows, columns).')
                # Angle bins: check the number of rows is equal to the number of angle bins minus one
                if weights.shape[0] != len(abins)-1:
                    raise ValueError(
                        'Number of weights rows must match number of angle bins minus one.')
                # Energy values: check the number of columns is equal to the number of energy values
                if weights.shape[1] != len(evalues):
                    raise ValueError(
                        'Number of weights columns must match number of energy values.')

                interpolation = angles['energy']['interpolation']
                # energy distribution type --> STILL TO HANDLE
                energy = openmc.stats.Tabular(
                    evalues, weights[i], interpolation=interpolation)
                # Handle strength
                strength = angles['strength']['data'][i]
                # Build source
                source = openmc.IndependentSource()
                source.paricle = particle
                source.space = space
                source.angle = angle
                source.energy = energy
                source.strength = strength
                # Append to source list
                angular_sources.append(source)
            # Append angular sources to the main source list
            sources.extend(angular_sources)

        return source

    def _build_tallies(self):
        
        tallies_data = self._benchmark_spec['tallies']

        # Initialize openmc tallies
        tallies = openmc.Tallies()

        for t in tallies_data:
            tally = openmc.Tally(name=t['name'])
            # Handle particle type
            particle = t['particle']
            particle_filter = openmc.ParticleFilter([particle])
            tally.filters.append(particle_filter)
            # Handle filters
            for d in t['filters']:
                if d['type'] == 'cell':
                    filter = openmc.CellFilter(d['values'])
                elif d['type'] == 'material':
                    filter = openmc.MaterialFilter(d['values'])
                elif d['type'] == 'surface':
                    filter = openmc.SurfaceFilter(d['values'])
                elif d['type'] == 'energy':
                    filter = openmc.EnergyFilter(d['values'])
                else:
                    raise ValueError(
                        f"Unsupported domain type: {d['type']}")

                tally.filters.append(filter)

            # Handle scores
            for q in t['scores']:
                tally.scores.append(q)

            # Store in tallies
            tallies.append(tally)

        return tallies

    def _build_settings(self):
        
        settings_data = self._benchmark_spec['settings']

        settings = openmc.Settings()
        if settings_data['run_mode'] == 'fixed_source':
            settings.run_mode = 'fixed source'
        elif settings_data['run_mode'] == 'k-eigenvalue':
            settings.run_mode = 'eigenvalue'
        else:
            raise ValueError(
                f"Unsupported run mode: {settings_data['run_mode']}")
        settings.batches = int(settings_data['batches'])
        settings.particles = int(settings_data['particles_per_batch'])
        settings.photon_transport = settings_data['photon_transport']
        # photon transport
        # weight windows
        # electron treatment
        settings.output = {'tallies': False}

        source = self._build_source()
        settings.source = source

        return settings

    def _build_model(self):
        
        materials = self._build_materials()
        geometry = self._build_geometry()
        settings = self._build_settings()
        tallies = self._build_tallies()
        model = openmc.Model(
            materials=materials,
            geometry=geometry,
            settings=settings,
            tallies=tallies
        )
        return model

    def _postprocess(self, statepoint: openmc.StatePoint | str | Path, mesh: str = 'mesh.h5m'):
        """Post-process the model after running."""
        # Retrieve tallies data from specifications
        tallies_data = self._benchmark_spec['tallies']
        
        tally_names = [t["name"] for t in tallies_data]
        normalizer = make_default_openmc_normalizer(mesh)

        # Accept both already-open StatePoint objects and statepoint file paths.
        code_version = "unknown"

        if isinstance(statepoint, openmc.StatePoint):
            sp = statepoint
            if hasattr(sp, "version"):
                code_version = ".".join(str(v) for v in sp.version)
            save_openmc_statepoint_tallies(
                statepoint=sp,
                filename="benchmark_results.h5",
                tally_names=tally_names,
                spec_tallies=tallies_data,
                tmc_coords={"realization": ["baseline"]},
                append_dim="realization",
                normalizer=normalizer,
            )
        else:
            with openmc.StatePoint(str(statepoint)) as sp:
                if hasattr(sp, "version"):
                    code_version = ".".join(str(v) for v in sp.version)
                save_openmc_statepoint_tallies(
                    statepoint=sp,
                    filename="benchmark_results.h5",
                    tally_names=tally_names,
                    spec_tallies=tallies_data,
                    tmc_coords={"realization": ["baseline"]},
                    append_dim="realization",
                    normalizer=normalizer,
                )

        if hasattr(self, "_write_spec_snapshot"):
            self._write_spec_snapshot("benchmark_results.h5")
        if hasattr(self, "_write_run_metadata"):
            self._write_run_metadata(
                code_name="openmc",
                code_version=code_version,
                nuclear_data_name=None,
                nuclear_data_version=None,
                geometry="cad",
                filename="benchmark_results.h5",
            )

        return

    def _uncertainty_quantification(self, *args, **kwargs):
        """Perform uncertainty quantification for the benchmark."""
        uq_data = self._benchmark_spec['uncertainty_quantification']
        tallies_data = self._benchmark_spec['tallies']

        mesh = 'mesh.h5m'

        # Run a TMC for every nuclide present in specifications
        for n, r in zip(uq_data[0]['nuclides'], uq_data[0]['realizations']):
            tmc_engine(model=self.model, realizations=r,
                       lib_name='endfb_80', nuclide=n[0], reaction=None,
                       perturb_xs=True, _is_benchmark=True, _spec_tallies=tallies_data,
                       _mesh=mesh, *args, **kwargs)

        return

    def run(
        self,
        uq: bool = False,
        *args,
        generate_report: bool = False,
        report_config: ReportConfig | None = None,
        **kwargs,
    ):
        """Run the benchmark simulation."""

        # Check if benchmark_results.h5 already exists and delete it
        path = Path("benchmark_results.h5")
        if path.exists():
            path.unlink()
            warnings.warn(
                f"Existing file '{path}' was found and deleted.", UserWarning)

        # Run the model
        if uq:
            # Perform uncertainty quantification
            self._uncertainty_quantification(*args, **kwargs)
        else:
            # Run the OpenMC model
            sp = self.model.run(*args, **kwargs)
            statepoint = openmc.StatePoint(sp)
            # Post-process the results
            self._postprocess(statepoint=statepoint)

        if generate_report:
            self._generate_report(report_config=report_config)

        return
