import networkx as nx
import pandas as pd

from D3mel.utilities.utilities_gos import UtilitiesGOs


def make_go_utilities():
    go = nx.MultiDiGraph()
    go.add_edge("GO:child", "GO:parent", key="is_a")
    go.add_edge("GO:parent", "GO:root", key="is_a")

    gaf_raw = pd.DataFrame(
        {
            "DB Object ID": ["FBgn0001"],
            "DB Object Symbol": ["gene1"],
            "GO ID": ["GO:child"],
            "DB Object Type": ["protein"],
            "Evidence Code": ["IMP"],
        }
    )
    return UtilitiesGOs(go, gaf_raw)


def test_get_ascendent_rejects_negative_depth():
    utilities = make_go_utilities()
    result = utilities.get_ascendent(["GO:child"], i_max=-1)
    assert result is None


def test_get_extend_handles_none_descendents():
    utilities = make_go_utilities()
    result = utilities.get_extend(["GO:child"], i_max_asc=1, i_max_desc=-1)
    assert result == ["GO:parent"]


def test_get_extend_handles_none_ascendents():
    utilities = make_go_utilities()
    result = utilities.get_extend(["GO:child"], i_max_asc=-1, i_max_desc=1)
    assert result == []
