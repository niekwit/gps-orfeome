import setuptools
import tomllib  # For Python >= 3.11
from pathlib import Path


# Function to read pyproject.toml
def read_pyproject_toml():
    pyproject_path = Path(__file__).parent / "pyproject.toml"
    if not pyproject_path.exists():
        raise FileNotFoundError(f"{pyproject_path} not found.")
    with open(pyproject_path, "rb") as f:
        return tomllib.load(f)


# Function to write versions to __init__.py
def write_version_info(package_name, version, container_image_version):
    init_file = Path(__file__).parent / "src" / package_name / "__init__.py"

    # Read existing content
    content = ""
    if init_file.exists():
        with open(init_file, "r") as f:
            content = f.read()

    # Check if __container_image_version__ already exists
    if "__container_image_version__" not in content:
        with open(init_file, "a") as f:
            f.write(f"\n__container_image_version__ = '{container_image_version}'\n")
    else:
        # Update existing line
        import re

        pattern = r'__container_image_version__\s*=\s*[\'"][^\'"]*[\'"]'
        replacement = f"__container_image_version__ = '{container_image_version}'"
        content = re.sub(pattern, replacement, content)
        with open(init_file, "w") as f:
            f.write(content)
    print(
        f"Injected __container_image_version__={container_image_version} into {init_file}"
    )


# Read data from pyproject.toml
pyproject_data = read_pyproject_toml()
package_name = pyproject_data["project"]["name"]
package_version = pyproject_data["project"]["version"]
container_image_version = (
    pyproject_data.get("tool", {})
    .get("gpsw", {})
    .get("container_image_version", "unknown")
)

# Write to __init__.py
write_version_info(package_name, package_version, container_image_version)

# Use setuptools.setup() as usual
setuptools.setup(
    # setuptools will now read project details from pyproject.toml
    # via build-backend = "setuptools.build_meta"
)
