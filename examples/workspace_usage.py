from flowpy import Workspace
import os

filename = "/home/swhite/Projects/flowpy_testing/test_data_2d_01.fcs"

wsp = Workspace("/home/swhite/Projects/flowpy_testing/20171218 Workspace.xml")

# First, get the gate IDs and names in the workspace.
# There are 2 possible classes of gates, one for each sample and then some
# for grouped samples
gates = wsp.get_gate_hierarchies()

# Apply a specified gate to a sample.
# Note, there are 2 API calls for analyzing samples, one for a sample gate
# and one for group gates
group_id = 2
results = wsp.analyze_sample(filename, None, 'group', group_id)

results['report'].to_csv(
    '%s__%s.csv' % (
        os.path.splitext(filename)[0],
        gates['groups'][str(group_id)]['name']),
    index=False,
    encoding='utf-8'
)
