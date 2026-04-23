from unittest.mock import MagicMock, patch

import pandas as pd

from D3mel.utilities.utilities import DatasetLoader, Utilities, _empty_optional_dataset


def test_empty_optional_dataset_gene_map_table():
    df = _empty_optional_dataset("gene_map_table")
    assert list(df.columns) == ["organism_abbreviation", "sequence_loc", "current_symbol"]


def test_empty_optional_dataset_isoforms():
    df = _empty_optional_dataset("dmel_unique_protein_isoforms")
    assert list(df.columns) == ["FB_gene_symbol"]


def test_dataset_loader_uses_public_fbd_api():
    with patch("D3mel.utilities.utilities.FBD") as mock_fbd:
        mock_instance = mock_fbd.return_value
        mock_instance.download_file.return_value = {"mock": "data"}

        loader = DatasetLoader()
        result = loader.load("go-basic")

        mock_instance.download_file.assert_called_once_with("go-basic")
        assert result == {"mock": "data"}


def test_dataset_loader_returns_empty_optional_dataset_on_failure():
    with patch("D3mel.utilities.utilities.FBD") as mock_fbd:
        mock_instance = mock_fbd.return_value
        mock_instance.download_file.side_effect = RuntimeError("download failed")

        loader = DatasetLoader()
        result = loader.load("gene_map_table", optional=True)

        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ["organism_abbreviation", "sequence_loc", "current_symbol"]


def test_dataset_loader_raises_for_required_dataset_on_failure():
    with patch("D3mel.utilities.utilities.FBD") as mock_fbd:
        mock_instance = mock_fbd.return_value
        mock_instance.download_file.side_effect = RuntimeError("download failed")

        loader = DatasetLoader()

        try:
            loader.load("go-basic", optional=False)
        except RuntimeError as exc:
            assert "download failed" in str(exc)
        else:
            raise AssertionError("Expected RuntimeError for required dataset failure")


def test_utilities_initialization_loads_expected_datasets():
    fake_go = MagicMock(name="go_graph")
    fake_gaf_raw = pd.DataFrame({"DB Object Symbol": [], "Evidence Code": [], "GO ID": [], "DB Object ID": [], "DB Object Type": []})
    fake_rpkm = pd.DataFrame({"Parent_library_name": [], "FBgn#": [], "RNASource_name": [], "RPKM_value": []})
    fake_gene_map = pd.DataFrame(columns=["organism_abbreviation", "sequence_loc", "current_symbol"])
    fake_isoforms = pd.DataFrame(columns=["FB_gene_symbol"])
    fake_gaf = pd.DataFrame({"DB Object ID": [], "DB Object Symbol": [], "GO ID": [], "DB Object Type": []})

    with patch.object(DatasetLoader, "load", side_effect=[fake_go, fake_gaf_raw, fake_rpkm, fake_gene_map, fake_isoforms]) as mock_load, \
         patch("D3mel.utilities.utilities.UtilitiesProteinComplex.__init__", return_value=None) as mock_pc_init, \
         patch("D3mel.utilities.utilities.UtilitiesGOs.__init__", return_value=None) as mock_gos_init, \
         patch("D3mel.utilities.utilities.UtilitiesPositionMap.__init__", return_value=None) as mock_pm_init, \
         patch("D3mel.utilities.utilities.UtilitiesRNAseq.__init__", return_value=None) as mock_rna_init, \
         patch("D3mel.utilities.utilities.UtilitiesGOs.get_GAF", return_value=fake_gaf) as mock_get_gaf:

        obj = Utilities()

        assert mock_load.call_count == 5
        mock_load.assert_any_call("go-basic")
        mock_load.assert_any_call("gene_association")
        mock_load.assert_any_call("gene_rpkm_report")
        mock_load.assert_any_call("gene_map_table", optional=True)
        mock_load.assert_any_call("dmel_unique_protein_isoforms", optional=True)

        mock_pc_init.assert_called_once_with(obj)
        mock_gos_init.assert_called_once_with(obj, fake_go, fake_gaf_raw)
        mock_pm_init.assert_called_once_with(obj, fake_gene_map, fake_isoforms)
        mock_rna_init.assert_called_once_with(obj, fake_rpkm)
        mock_get_gaf.assert_called_once_with(obj)
        assert obj.gaf is fake_gaf
