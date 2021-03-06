from xml.etree import cElementTree
import os
from flowpy import Sample
from flowpy.models import gate
from copy import deepcopy


def parse_gate_element(gate_element):
    gate_type, _ = gate_element.tag.split('Gate')

    if gate_type == 'Polygon':
        params = []

        params_el = gate_element.find('ParameterNames').find('StringArray').findall('String')

        for param in params_el:
            params.append(param.text)

        # PolygonGate can contain either a Polygon or PolyRect element, and the attributes
        # of this element reveal the parameter used for the x & y axes
        if len(gate_element.findall('Polygon')) > 0:
            poly_el = gate_element.find('Polygon')
        elif len(gate_element.findall('PolyRect')) > 0:
            poly_el = gate_element.find('PolyRect')
        else:
            return {}

        x_axis = poly_el.attrib['xAxisName']
        y_axis = poly_el.attrib['yAxisName']

        vertex_els = gate_element.iter('Vertex')
        vertices = []

        for v in vertex_els:
            vertices.append(v.attrib)

        gate_dict = {
            'type': gate_type,
            'parameters': params,
            'x_axis': x_axis,
            'y_axis': y_axis,
            'vertices': vertices
        }
    elif gate_type == 'Boolean':
        groups = []

        gate_paths = gate_element.find('GatePaths').iter('String')

        for g in gate_paths:
            groups.append(g.text)

        gate_dict = {
            'type': gate_type,
            'specification': gate_element.attrib['specification'],
            'groups': groups
        }
    else:
        gate_dict = {}

    return gate_dict


def find_group_samples(group_element):
    group_samples = group_element.iter('SampleRef')

    samples = []

    for gs in group_samples:
        samples.append(gs.attrib['sampleID'])

    return samples


def find_nested_populations(element):
    populations = element.findall('Population')

    population_dict = {}

    for pop in populations:
        gates = []

        for child in pop:
            if 'Gate' in child.tag:
                gates.append(parse_gate_element(child))

        child_pops = find_nested_populations(pop)

        population_dict[pop.attrib['nodeName']] = {
            'gates': gates,
            'children': child_pops
        }

    return population_dict


class Workspace(object):
    """
    Parses an FlowJo XML workspace file to extract gating information
    """
    def __init__(self, xml_workspace):
        """
        
        :param xml_workspace: FlowJo XML workspace file
        """

        with open(xml_workspace, 'r') as in_file:
            tree = cElementTree.parse(in_file)

        # first parse out the samples
        samples = tree.iter('Sample')

        sample_dict = {}

        for s in samples:
            sample_node = s.find('SampleNode')
            filename = sample_node.attrib['nodeName']

            if filename == "":
                continue

            # a sample may have it's own individual gating hierarchy independent of
            # a sample group, we'll search for them here
            sample_populations = find_nested_populations(sample_node)

            sample_dict[s.attrib['sampleID']] = {
                'filename': sample_node.attrib['nodeName'],
                'eventCount': s.attrib['eventCount'],
                'populations': sample_populations
            }

        # next, parse the sample groups
        groups = tree.iter('GroupNode')

        group_dict = {}

        for g in groups:
            pops = find_nested_populations(g)

            # don't store groups with no gates
            if len(pops) <= 0:
                continue

            group_samples = find_group_samples(g)

            group_dict[g.attrib['groupID']] = {
                'group_name': g.attrib['nodeName'],
                'samples': group_samples,
                'populations': pops
            }

        # TODO: add comp attribute for lookup
        self.samples = sample_dict
        self.groups = group_dict

    def get_gate_hierarchies(self):
        gate_dict = {
            'samples': {},
            'groups': {}
        }

        for s_id, s_dict in self.samples.items():
            gate_dict['samples'][s_id] = {
                'name': s_dict['filename']
            }

        for g_id, g_dict in self.groups.items():
            gate_dict['groups'][g_id] = {
                'name': g_dict['group_name'],
                'samples': g_dict['samples']
            }

        return gate_dict

    def find_matching_gate_hierarchies(self, fcs_file_path):
        base_name = os.path.basename(fcs_file_path)
        chosen_sample = None

        for sample_id, ws_sample in self.samples.items():
            if base_name == ws_sample['filename']:
                chosen_sample = sample_id
                break

        if chosen_sample is None:
            UserWarning("%s was not found in workspace" % base_name)
            return None

        matching_gates = {
            'sample_id': chosen_sample,
            'groups': []
        }

        # find groups to which the sample belongs
        for group_id, ws_group in self.groups.items():
            if chosen_sample in ws_group['samples']:
                matching_gates['groups'].append(group_id)

        return matching_gates

    def analyze_sample(self, fcs_file_path, comp_matrix, gate_type, gate_id):
        base_name = os.path.basename(fcs_file_path)

        s = Sample(fcs_file_path, track_indices=True)
        s.generate_subsample(0, random_seed=123)
        s.compensate_events(comp_matrix)
        
        # check gate type, options are 'sample' or 'group'
        if gate_type == 'sample':
            chosen_gate = deepcopy(self.samples[str(gate_id)])
        elif gate_type == 'group':
            chosen_gate = deepcopy(self.groups[str(gate_id)])
        else:
            raise ValueError("Gate type %s is not valid, use 'sample' or 'group'" % gate_type)

        results_dict = {
            'filename': base_name
        }

        # looks like the FlowJo XML gates are saved
        # on compensated but not transformed data
        gate.apply_gating_hierarchy(
            s.events_compensated,
            s.channels,
            chosen_gate['populations']
        )

        results_dict['populations'] = chosen_gate['populations']
        results_dict['report'] = gate.results_to_dataframe(chosen_gate)

        return results_dict
