from flowpy import Workspace
import pandas as pd
import os

#wsp = Workspace("example_workspace.xml")
wsp = Workspace("/home/swhite/Projects/HMMVB/data/EQAPOL_normal/final_comparison/flowjo_workspace/Denny ICS 06Nov09lym_reexp.xml")

comp_matrix_1 = """Blue B-A	Blue A-A	Violet H-A	Violet G-A	Red C-A	Red B-A	Red A-A	Green E-A	Green D-A	Green C-A	Green A-A
1	7.5e-3	0	0.0414	0	0	0	9.e-3	7.5e-3	0	0
-0.01	1	0	0.035	0.02	0.12	-0.045	0	0.01	0.145	0.055
-0.02	-0.01	1	1.1	0	0	0	0	0	0	0.05
0.185	5.17e-4	0.0229	1	0	-0.01	0	0	0	-0.04	0
8.e-3	0	-5.e-3	0.1	1	0.471	0.15	-0.01	1.08e-3	0.4	0.0229
0	0.0517	1.e-3	0	0	1	0.579	0	3.7e-3	0.03	0.058
0	0	0	0	0.0258	0	1	0	0	0.0107	0.136
0	0.0364	0	-5.e-3	0	0	0	1	0.808	0.152	9.e-3
5.e-3	0.06	-1.e-3	0	-5.e-3	5.e-3	-0.01	0.0584	1	0.265	0.015
0	0.29	0	0	0.091	0.084	0.028	0.0162	0.0174	1	0.093
0	0.01	-5.e-4	-5.e-3	-1.e-3	0	0.0885	0.023	0.0225	4.53e-3	1"""

comp_matrix_2 = """Blue B-A	Blue A-A	Violet H-A	Violet G-A	Red C-A	Red B-A	Red A-A	Green E-A	Green D-A	Green C-A	Green A-A
1	7.5e-3	0	0.0414	0	0	0	9.e-3	7.5e-3	0	0
-0.01	1	0	0.035	0.02	0.12	-0.045	0	0.01	0.145	0.055
-0.02	-0.25	1	1.1	0	0	0	0	0	0	0.04
0.185	0.02	0.0229	1	0	-0.01	0	0	0	-0.04	0
8.e-3	0	-5.e-3	0.1	1	0.471	0.15	-0.01	1.08e-3	0.4	0.0229
0	0.0517	1.e-3	0	0	1	0.579	0	3.7e-3	0.03	0.058
0	0	0	0	0.0258	0	1	0	0	0.0107	0.136
0	0.0364	0	-5.e-3	0	0	0	1	0.808	0.152	9.e-3
5.e-3	0.06	-1.e-3	-0.01	-5.e-3	-2.5e-3	-0.01	0.0584	1	0.265	0.015
0	0.29	0	0	0.091	0.08	0.028	0.0162	0.0174	1	0.093
0	0.01	-5.e-4	-5.e-3	-1.e-3	0	0.0885	0.023	0.0225	4.53e-3	1"""

filename = "F69018CN-01_SEB.fcs"

results = wsp.analyze_sample(filename, comp_matrix_1)


def parse_results_dict(population_root, parent_label):
    parsed_results = []

    for label, pop in population_root.items():
        if len(pop['gates']) > 1:
            print('multi-region gate')

        parsed_results.append({
            'parents': parent_label,
            'label': label,
            'type': pop['gates'][0]['type'],
            'count': pop['gates'][0]['result']['gated_events'].shape[0],
            'parent_count': pop['gates'][0]['result']['ungated_count']
        })

        child_results = parse_results_dict(pop['children'], ' -> '.join([parent_label, label]))
        parsed_results.extend(child_results)

    return parsed_results


def results_to_dataframe(results):
    for sg in results['sample_groups']:
        sg_results = parse_results_dict(sg['populations'], 'root')
        sg_df = pd.DataFrame(sg_results)
        sg_df.sort(columns=['parents', 'label'], inplace=True)
        sg_df.to_csv(
            '%s__%s.csv' % (os.path.splitext(filename)[0], sg['group_name']),
            columns=['parents', 'label', 'type', 'parent_count', 'count'],
            index=False,
            encoding='utf-8'
        )

results_to_dataframe(results)
