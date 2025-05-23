import unittest
from unittest.mock import patch, MagicMock, mock_open
import subprocess
from pathlib import Path


# Import the functions from your utils module
# Assuming your utils.py is in src/gpsw/utils.py
# and your test file is at the project root or in a 'tests' directory
# You might need to adjust this import based on your actual project structure
# If running tests from the root, ensure the package is installed -e or use python -m unittest
from src.gpsw import utils


class TestUtils(unittest.TestCase):

    # --- Test get_package_version ---
    @patch("src.gpsw.utils.version")
    def test_get_package_version_installed(self, mock_version):
        mock_version.return_value = (
            "1.2.3"  # This is the mocked return value for this specific test
        )
        self.assertEqual(utils.get_package_version(), "1.2.3")
        mock_version.assert_called_once_with("gpsw")

    # Correct patch target for 'version' and 'PackageNotFoundError'
    @patch(
        "src.gpsw.utils.version", side_effect=utils.PackageNotFoundError
    )  # Use utils.PackageNotFoundError
    def test_get_package_version_not_installed(self, mock_version):
        self.assertEqual(utils.get_package_version(), "v0.0.0-dev")
        mock_version.assert_called_once_with("gpsw")

    # --- Test get_latest_release_tag ---
    @patch("src.gpsw.utils.Github")  # Patch the Github import within your utils module
    def test_get_latest_release_tag_success(self, MockGithub):
        mock_repo = MagicMock()
        mock_release = MagicMock()
        mock_release.tag_name = "v1.0.0"
        mock_repo.get_latest_release.return_value = mock_release
        MockGithub.return_value.get_user.return_value.get_repo.return_value = mock_repo

        self.assertEqual(utils.get_latest_release_tag(), "v1.0.0")

    @patch("src.gpsw.utils.Github")
    @patch("builtins.print")  # Mock print to check output
    def test_get_latest_release_tag_failure(self, mock_print, MockGithub):
        MockGithub.side_effect = Exception("Network error")  # Simulate a network error
        self.assertIsNone(utils.get_latest_release_tag())
        mock_print.assert_called_with(
            "Error fetching release using PyGithub: Network error"
        )

    # --- Test dry_run ---
    @patch("subprocess.run")
    @patch("builtins.print")
    @patch("sys.exit")
    def test_dry_run_success_quiet_false(
        self, mock_exit, mock_print, mock_subprocess_run
    ):
        # Arrange
        mock_args = MagicMock()
        mock_args.quiet = False
        mock_subprocess_run.return_value = MagicMock(returncode=0)

        # Act
        utils.dry_run(mock_args)

        # Assert
        mock_subprocess_run.assert_called_once_with(
            ["snakemake", "-n", "-p"], check=True
        )
        mock_print.assert_any_call("Performing dry-run of the workflow")
        mock_print.assert_any_call("Dry-run completed successfully!")
        mock_exit.assert_called_once_with(0)

    @patch("subprocess.run")
    @patch("builtins.print")
    @patch("sys.exit")
    def test_dry_run_success_quiet_true(
        self, mock_exit, mock_print, mock_subprocess_run
    ):
        # Arrange
        mock_args = MagicMock()
        mock_args.quiet = True
        mock_subprocess_run.return_value = MagicMock(returncode=0)

        # Act
        utils.dry_run(mock_args)

        # Assert
        mock_subprocess_run.assert_called_once_with(
            ["snakemake", "-n", "--quiet", "all"], check=True
        )
        mock_print.assert_any_call("Performing dry-run of the workflow")
        mock_print.assert_any_call("Dry-run completed successfully!")
        mock_exit.assert_called_once_with(0)

    @patch("subprocess.run", side_effect=subprocess.CalledProcessError(1, "cmd"))
    @patch("builtins.print")
    @patch("sys.exit")
    def test_dry_run_failure_quiet_false(
        self, mock_exit, mock_print, mock_subprocess_run
    ):
        # Arrange
        mock_args = MagicMock()
        mock_args.quiet = False

        # Act
        utils.dry_run(mock_args)

        # Assert
        mock_subprocess_run.assert_called_once_with(
            ["snakemake", "-n", "-p"], check=True
        )
        mock_print.assert_any_call("Performing dry-run of the workflow")
        mock_exit.assert_called_once_with(1)

    @patch("subprocess.run")
    @patch("builtins.print")
    @patch("sys.exit")
    def test_dry_run_failure_quiet_true_reruns_verbose(
        self, mock_exit, mock_print, mock_subprocess_run
    ):
        # Arrange
        mock_args = MagicMock()
        mock_args.quiet = True

        # Set side_effect to be a list:
        # The first call to subprocess.run will raise CalledProcessError
        # The second call will return a successful mock object (like a normal completion)
        mock_subprocess_run.side_effect = [
            subprocess.CalledProcessError(
                1, "snakemake -n --quiet all"
            ),  # First call fails
            MagicMock(returncode=0),  # Second call succeeds
        ]

        # Act
        utils.dry_run(mock_args)

        # Assert
        # Check that subprocess.run was called twice
        self.assertEqual(mock_subprocess_run.call_count, 2)

        # The first call (failed quiet dry-run)
        mock_subprocess_run.call_args_list[0].assert_called_with(
            ["snakemake", "-n", "--quiet", "all"], check=True
        )
        # The second call (verbose dry-run for error message)
        mock_subprocess_run.call_args_list[1].assert_called_with(["snakemake", "-np"])

        mock_print.assert_any_call("Performing dry-run of the workflow")
        # No 'Dry-run completed successfully!' should be printed for a failure case
        # You might want to check the absence of this print too if desired
        mock_exit.assert_called_once_with(1)

    # --- Test create_rule_graph ---
    @patch("os.makedirs")
    @patch("subprocess.run")
    @patch("pydot.graph_from_dot_file")
    def test_create_rule_graph_success(
        self, mock_graph_from_dot_file, mock_subprocess_run, mock_makedirs
    ):
        # Arrange
        mock_graph = MagicMock()
        mock_graph_from_dot_file.return_value = [mock_graph]  # pydot returns a list

        # Act
        utils.create_rule_graph()

        # Assert
        mock_makedirs.assert_called_once_with("images", exist_ok=True)
        expected_command_part = r"snakemake --quiet all --forceall --rulegraph | grep -v '\-> 0\|0\[label = \"all\"' > images/rulegraph.dot"
        mock_subprocess_run.assert_called_once_with(expected_command_part, shell=True)
        mock_graph_from_dot_file.assert_called_once_with("images/rulegraph.dot")
        mock_graph.write_pdf.assert_called_once_with("images/rulegraph.pdf")

    # --- Test profile_arg ---
    # Mock os.path.exists and configparser behavior
    @patch("os.path.exists")
    @patch("configparser.ConfigParser")
    @patch("builtins.open", new_callable=mock_open)  # Mock file open
    @patch("os.makedirs")  # Mock makedirs for safety, though not called here
    @patch("pathlib.Path.home")  # <--- IMPORTANT: Mock Path.home()
    def test_profile_arg_none_profile(
        self,
        mock_path_home,
        MockOsMakedirs,
        mock_file_open,
        MockConfigParser,
        mock_exists,
    ):
        # Configure mock_path_home to return a predictable home directory
        mock_path_home.return_value = Path("/mock/home/niek")

        mock_args = MagicMock(profile="None")
        self.assertEqual(utils.profile_arg(mock_args), [""])
        mock_exists.assert_not_called()  # No file operations if profile is "None"

    @patch("os.path.exists", return_value=True)  # Config file exists
    @patch("configparser.ConfigParser")
    @patch("builtins.open", new_callable=mock_open)  # Mock file open
    @patch("os.makedirs")  # Mock makedirs for safety, though not called here
    @patch("pathlib.Path.home")  # <--- IMPORTANT: Mock Path.home()
    def test_profile_arg_config_exists_profile_in_config(
        self,
        mock_path_home,
        MockOsMakedirs,
        mock_file_open,
        MockConfigParser,
        mock_exists,
    ):
        # Configure mock_path_home to return a predictable home directory
        mock_path_home.return_value = Path("/mock/home/niek")
        mock_config_file_path = Path("/mock/home/niek/.gpsw/config.ini")

        mock_args = MagicMock(profile=None)  # User didn't provide it
        mock_config = MagicMock()
        mock_config.__contains__.side_effect = (
            lambda key: key == "DEFAULT"
        )  # Mock 'DEFAULT' section check
        mock_config.__getitem__.return_value = {
            "profile": "/path/from/config"
        }  # Mock profile value
        MockConfigParser.return_value = mock_config

        # Act
        result = utils.profile_arg(mock_args)

        # Assert
        self.assertEqual(result, ["--profile", "/path/from/config"])
        mock_exists.assert_called_once_with(
            mock_config_file_path
        )  # <--- Use the Path object here
        mock_config.read.assert_called_once_with(
            mock_config_file_path
        )  # <--- Use the Path object here
        mock_file_open.assert_not_called()  # No file written
        MockOsMakedirs.assert_not_called()

    @patch("os.path.exists", return_value=True)  # Config file exists
    @patch("configparser.ConfigParser")
    @patch("builtins.open", new_callable=mock_open)  # Mock file open for writing
    @patch("os.makedirs")  # Mock makedirs for safety, though not called here
    @patch("pathlib.Path.home")  # <--- IMPORTANT: Mock Path.home()
    def test_profile_arg_config_exists_profile_not_in_config_but_in_args(
        self,
        mock_path_home,
        MockOsMakedirs,
        mock_file_open,
        MockConfigParser,
        mock_exists,
    ):
        # Configure mock_path_home to return a predictable home directory
        mock_path_home.return_value = Path("/mock/home/niek")
        mock_config_file_path = Path("/mock/home/niek/.gpsw/config.ini")

        mock_args = MagicMock(profile="/path/from/args")  # User provided it
        mock_config = MagicMock()
        mock_config.__contains__.side_effect = lambda key: key == "DEFAULT"
        mock_config.__getitem__.side_effect = lambda key: (
            {} if key == "DEFAULT" else None
        )  # No 'profile' in 'DEFAULT'
        MockConfigParser.return_value = mock_config

        # Act
        result = utils.profile_arg(mock_args)

        # Assert
        self.assertEqual(result, ["--profile", "/path/from/args"])
        mock_exists.assert_called_once_with(
            mock_config_file_path
        )  # <--- Use the Path object here
        mock_config.read.assert_called_once_with(
            mock_config_file_path
        )  # <--- Use the Path object here
        mock_config.__setitem__.assert_called_once_with(
            "DEFAULT", {"profile": "/path/from/args"}
        )
        mock_file_open.assert_called_once_with(
            mock_config_file_path, "w"
        )  # <--- Use the Path object here
        mock_config.write.assert_called_once_with(mock_file_open())
        MockOsMakedirs.assert_not_called()

    @patch("os.path.exists", return_value=False)  # Config file does NOT exist
    @patch("configparser.ConfigParser")
    @patch("builtins.open", new_callable=mock_open)  # Mock file open for writing
    @patch("os.makedirs")  # Mock makedirs for creating ~/.gpsw
    @patch("pathlib.Path.home")  # <--- IMPORTANT: Mock Path.home()
    def test_profile_arg_config_does_not_exist_profile_in_args(
        self,
        mock_path_home,
        mock_makedirs,
        mock_file_open,
        MockConfigParser,
        mock_exists,
    ):
        # Configure mock_path_home to return a predictable home directory
        mock_path_home.return_value = Path("/mock/home/niek")
        mock_config_dir = Path("/mock/home/niek/.gpsw")
        mock_config_file_path = Path("/mock/home/niek/.gpsw/config.ini")

        mock_args = MagicMock(profile="/path/from/new_args")
        mock_config = MagicMock()
        MockConfigParser.return_value = mock_config

        # Act
        result = utils.profile_arg(mock_args)

        # Assert
        self.assertEqual(result, ["--profile", "/path/from/new_args"])
        mock_exists.assert_called_once_with(
            mock_config_file_path
        )  # <--- Use the Path object here
        mock_makedirs.assert_called_once_with(
            mock_config_dir, exist_ok=True
        )  # <--- Use the Path object here
        mock_config.__setitem__.assert_called_once_with(
            "DEFAULT", {"profile": "/path/from/new_args"}
        )
        mock_file_open.assert_called_once_with(
            mock_config_file_path, "w"
        )  # <--- Use the Path object here
        mock_config.write.assert_called_once_with(mock_file_open())

    @patch("os.path.exists", return_value=False)  # Config file does NOT exist
    @patch("configparser.ConfigParser")
    @patch("builtins.open", new_callable=mock_open)
    @patch("os.makedirs")
    @patch("pathlib.Path.home")  # <--- IMPORTANT: Mock Path.home()
    def test_profile_arg_config_does_not_exist_no_profile_in_args(
        self,
        mock_path_home,
        MockOsMakedirs,
        mock_file_open,
        MockConfigParser,
        mock_exists,
    ):
        # Configure mock_path_home to return a predictable home directory
        mock_path_home.return_value = Path("/mock/home/niek")
        mock_config_file_path = Path("/mock/home/niek/.gpsw/config.ini")

        mock_args = MagicMock(profile=None)

        # Act & Assert
        with self.assertRaises(ValueError) as cm:
            utils.profile_arg(mock_args)
        self.assertIn("No profile provided at first run of GPSW.", str(cm.exception))

        mock_exists.assert_called_once_with(
            mock_config_file_path
        )  # <--- Use the Path object here
        MockOsMakedirs.assert_not_called()
        mock_file_open.assert_not_called()


if __name__ == "__main__":
    unittest.main()
