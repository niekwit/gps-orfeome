import unittest

# import os
# import tempfile
# import shutil
from unittest.mock import patch, MagicMock, mock_open
import subprocess
from pathlib import Path
import builtins

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
        # Set captured_output and text to False for non-quiet run to match function's default behavior
        mock_subprocess_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        # Act
        utils.dry_run(mock_args)

        # Assert
        # capture_output is False when not quiet (default stream behavior)
        mock_subprocess_run.assert_called_once_with(
            ["snakemake", "-n", "-p"], check=True, capture_output=False, text=True
        )
        mock_print.assert_any_call("Performing dry-run of the workflow")
        mock_print.assert_any_call(
            "Dry-run completed successfully!"
        )  # This print is now unconditional on success
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
        # Ensure stdout and stderr are captured and empty for a successful quiet run
        mock_subprocess_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        # Act
        utils.dry_run(mock_args)

        # Assert
        # Expecting 'all' argument with --quiet as per original function logic
        mock_subprocess_run.assert_called_once_with(
            ["snakemake", "-n", "--quiet", "all"],
            check=True,
            capture_output=True,
            text=True,
        )
        mock_print.assert_any_call("Performing dry-run of the workflow")
        mock_print.assert_any_call(
            "Dry-run completed successfully!"
        )  # This print is now unconditional on success
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
            ["snakemake", "-n", "-p"], check=True, capture_output=False, text=True
        )
        mock_print.assert_any_call("Performing dry-run of the workflow")
        mock_print.assert_any_call(
            unittest.mock.ANY
        )  # Check for any error message print
        mock_exit.assert_called_once_with(1)

    @patch("builtins.print")
    @patch("subprocess.run")
    def test_fn(self, mock_subprocess_run, mock_print):
        mock_args = MagicMock()
        mock_args.quiet = True

        mock_called_process_error = subprocess.CalledProcessError(
            1,
            ["snakemake", "-n", "--quiet", "all"],
            stderr="ERROR: Something went wrong in Snakemake\n",
        )

        # First call fails, second is verbose fallback
        mock_subprocess_run.side_effect = [
            mock_called_process_error,
            MagicMock(returncode=0, stdout="Verbose run output", stderr=""),
        ]

        with self.assertRaises(SystemExit) as cm:
            utils.dry_run(mock_args)

        self.assertEqual(cm.exception.code, 1)

        self.assertEqual(mock_subprocess_run.call_count, 2)

        printed_lines = [call.args[0] for call in mock_print.call_args_list]

        print("\nCaptured prints:")
        for line in printed_lines:
            print(f"LINE: {repr(line)}")

        assert any(
            ("Attempting re-run" in line and "--quiet" in line)
            for line in printed_lines
        ), "Expected a retry message with '--quiet' in the output"

    # --- Test create_rule_graph ---
    # This test needed significant refactoring due to how category_colors is a local variable.
    # The fix is to calculate expected_category_colors within the test itself, mirroring the function's logic.
    @patch("src.gpsw.utils.extract_rule_categories")
    @patch("pydot.graph_from_dot_file")
    @patch("builtins.open", new_callable=mock_open)
    @patch("subprocess.run")
    @patch("os.makedirs")
    def test_create_rule_graph_success(
        self,
        mock_makedirs,
        mock_subprocess_run,
        mock_file_open,
        mock_graph_from_dot_file,
        mock_extract_categories,
    ):
        # Arrange
        # Mock extract_rule_categories return value with a few categories for comprehensive testing
        mock_extract_categories.return_value = {
            "cutadapt": "Preprocessing",
            "fastqc": "Quality Control",
            "align": "Alignment",
            "merge": "Post-processing",
        }

        # DOT output from snakemake. Include 'all' rule (node 0) for filtering test
        raw_snakemake_stdout = """
digraph G {
    0[label = "all"];
    1[label = "cutadapt"];
    2[label = "fastqc"];
    3[label = "align"];
    4[label = "merge"];
    2 -> 1;
    1 -> 0; // Edge to 'all' - should be filtered
    3 -> 1;
    4 -> 3;
}
"""
        mock_subprocess_run.return_value = MagicMock(stdout=raw_snakemake_stdout)

        mock_graph = MagicMock()
        mock_graph_from_dot_file.return_value = [mock_graph]

        # Act
        utils.create_rule_graph()

        # Assert
        mock_makedirs.assert_called_once_with("images", exist_ok=True)

        mock_subprocess_run.assert_called_once_with(
            ["snakemake", "--quiet", "all", "--forceall", "--rulegraph"],
            capture_output=True,
            text=True,
            check=True,
        )

        mock_file_open.assert_any_call("images/rulegraph.dot", "w")
        handle = mock_file_open()
        written = "".join(call.args[0] for call in handle.write.call_args_list)

        # Manually determine the expected colors based on the function's logic
        palette = [
            "#F0FAF0",
            "#FAF0F0",
            "#F0F0FA",
            "#FFF8E1",
            "#E6F0FF",
            "#FDEDEC",
            "#F3F3F3",
        ]

        # This mirrors the logic inside create_rule_graph for category_colors
        # Ensure categories are sorted alphabetically to match the function's behavior
        # Get unique categories from the mocked extract_rule_categories return value
        unique_categories_from_mock = list(
            set(mock_extract_categories.return_value.values())
        )
        sorted_unique_categories = sorted(unique_categories_from_mock)

        expected_category_colors = {
            category: palette[i % len(palette)]
            for i, category in enumerate(sorted_unique_categories)
        }

        # Assert on key parts of the written content to ensure correct structure and colors
        self.assertIn("digraph snakemake_dag {", written)
        self.assertIn("rankdir=TB;", written)
        self.assertIn(
            'node [shape=box, style="rounded,filled", fontname=Helvetica, fontsize=10, penwidth=1.5, color=black, fillcolor=white];',
            written,
        )
        self.assertIn('edge [color="#888888", penwidth=1.2];', written)

        # Check for each subgraph cluster by iterating through the expected categories
        # and checking for both label and fillcolor.
        # Note: Order of subgraphs in the written output depends on sorted(clusters) in the function.
        # This will be based on alphabetical order of categories.

        # Expected categories (alphabetical order): 'Alignment', 'Post-processing', 'Preprocessing', 'Quality Control'

        # Alignment
        self.assertIn(f"subgraph cluster_alignment {{", written)
        self.assertIn(f'label = "Alignment";', written)
        self.assertIn(
            f'fillcolor = "{expected_category_colors["Alignment"]}";', written
        )
        self.assertIn(f'3 [label="align"];', written)

        # Post-processing
        self.assertIn(f"subgraph cluster_post_processing {{", written)
        self.assertIn(f'label = "Post-processing";', written)
        self.assertIn(
            f'fillcolor = "{expected_category_colors["Post-processing"]}";', written
        )
        self.assertIn(f'4 [label="merge"];', written)

        # Preprocessing
        self.assertIn(f"subgraph cluster_preprocessing {{", written)
        self.assertIn(f'label = "Preprocessing";', written)
        self.assertIn(
            f'fillcolor = "{expected_category_colors["Preprocessing"]}";', written
        )
        self.assertIn(f'1 [label="cutadapt"];', written)

        # Quality Control
        self.assertIn(f"subgraph cluster_quality_control {{", written)
        self.assertIn(f'label = "Quality Control";', written)
        self.assertIn(
            f'fillcolor = "{expected_category_colors["Quality Control"]}";', written
        )
        self.assertIn(f'2 [label="fastqc"];', written)

        # Check edges, ensuring 'all' edge is filtered
        self.assertIn("2 -> 1;", written)
        self.assertIn("3 -> 1;", written)
        self.assertIn("4 -> 3;", written)
        self.assertNotIn("-> 0;", written)  # Crucially, 'all' edge should be gone
        self.assertNotIn('0[label = "all"', written)  # 'all' node should be gone

        self.assertIn("}", written)  # Ensure closing bracket is there

        mock_graph_from_dot_file.assert_called_once_with("images/rulegraph.dot")
        mock_graph.write_pdf.assert_called_once_with("images/rulegraph.pdf")
        mock_extract_categories.assert_called_once_with("workflow/rules")

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
