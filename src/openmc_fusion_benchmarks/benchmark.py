import yaml
from pathlib import Path
from abc import ABC, abstractmethod
import openmc
import numpy as np
from .validate import validate_benchmark


BENCHMARK_DIR = Path(__file__).parent / "benchmarks"


class Benchmark(ABC):
    def __init__(self, name: str):
        self.name = name

        # Validate the benchmark specification
        validate_benchmark(name)

        with (BENCHMARK_DIR / f"{self.name}/specifications.yaml").open("r") as f:
            self._benchmark_spec = yaml.safe_load(f)

        self._read_metadata()

    @abstractmethod
    def build_materials(self):
        """Build materials for the benchmark."""
        pass

    @abstractmethod
    def build_geometry(self):
        """Build geometry for the benchmark."""
        pass

    @abstractmethod
    def build_source(self):
        """Build source for the benchmark."""
        pass

    @abstractmethod
    def build_settings(self):
        """Build settings for the benchmark."""
        pass

    @abstractmethod
    def build_tallies(self):
        """Build tallies for the benchmark."""
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


class OpenmcBenchmark(Benchmark):
    def __init__(self, name: str):
        super().__init__(name)
        self._materials = None
        self._geometry = None
        self._settings = None
        self._tallies = None

    def build_materials(self):
        # Implement the logic to build materials for OpenMC
        material_data = self._benchmark_spec['materials']

        fraction_map = {'atomic': 'ao', 'weight': 'wo'}

        materials = openmc.Materials()
        for m in material_data:
            mat = openmc.Material(name=m['name'])
            mat.id = m['id']
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

    def build_geometry(self):
        raise NotImplementedError(
            'Geometry building not implemented yet. Please implement the build_geometry method in the OpenmcBenchmark class.'
        )

    def build_source(self):
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
            abins = np.array(angles['angle']['bins'])
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

    def build_tallies(self):
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

    def build_settings(self):
        settings_data = self._benchmark_spec['settings']

        settings = openmc.Settings()
        if settings_data['run_mode'] == 'fixed source':
            settings.run_mode = 'fixed source'
        elif settings_data['run_mode'] == 'k-eigenvalue':
            settings.run_mode = 'eigenvalue'
        settings.batches = int(settings_data['batches'])
        settings.particles = int(settings_data['particles_per_batch'])
        settings.photon_transport = settings_data['photon_transport']
        # photon transport
        # weight windows
        # electron treatment
        settings.output = {'tallies': False}

        source = self.build_source()
        settings.source = source

        return settings
