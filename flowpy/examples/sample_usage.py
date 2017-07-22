from flowpy import Sample

# import numpy as np
#
# a = np.array(
#     [
#         [0, 1, 2, 3],
#         [0, 1, 2, 3],
#         [0, 1, 2, 3],
#         [0, 1, 2, 3],
#         [0, 1, 2, 3],
#     ],
#     dtype=np.float
# )
#
# indices = [3, 2]
# new_a = a.copy()[:, indices]
# new_a += 0.1
# a[:, indices] = new_a
# print(new_a.shape)

comp_matrix = """Blue B-A	Blue A-A	Violet H-A	Violet G-A	Red C-A	Red B-A	Red A-A	Green E-A	Green D-A	Green C-A	Green A-A
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

s = Sample("K690194C-01_Costim.fcs")
s.generate_subsample(1000, random_seed=123)
s.compensate_events(comp_matrix)
s.transform_logicle()

print(s)

