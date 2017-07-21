from matplotlib.path import Path
import numpy as np
import warnings


def parse_polygon_gate(events, channel_labels, gate):
    """
    Extract events in given Polygon gate

    :param events: NumPy array of events on which to apply the gate
    :param channel_labels: dictionary of channel labels (keys are channel #'s)
    :param gate: dictionary for a 'Polygon' gate
    :return:
    """
    # First, get the column indices for the x and y parameters
    x_label = gate['x_axis']
    y_label = gate['y_axis']
    x_index = None
    y_index = None

    for chan_number, labels in channel_labels.items():
        if labels['PnN'] in x_label:
            x_index = int(chan_number) - 1
        elif labels['PnN'] in y_label:
            y_index = int(chan_number) - 1

    if x_index is None or y_index is None:
        raise ValueError("Channel labels not found in data for polygon gate")

    xy_events = events[:, [x_index, y_index]]

    xy_vertices = []

    for vertex in gate['vertices']:
        xy_vertices.append(
            [
                float(vertex['x']),
                float(vertex['y'])
            ]
        )

    path = Path(xy_vertices)

    is_in_gate = path.contains_points(xy_events, radius=0.0)

    gated_events = events[is_in_gate]

    return {
        'gated_events': gated_events,
        'ungated_count': events.shape[0]
    }


def parse_gates(events, channel_labels, gating_dict):
    # We ignore boolean gates here since they reference other gates that need to be
    # calculated first. The boolean gates will be handled by another function that
    # should be run after this one. Don't really like this as it means at least 2
    # passes through the gating_dict, but it will work for now.
    boolean_gates = []

    for gate_label, gate_objects in gating_dict.items():
        # Iterate through gates, I suppose it is possible to have multiple sub-regions
        # drawn for a single "gate" so we'll iterate and check in case there's > 1 gate
        gated_events = np.ndarray(shape=(0, events.shape[1]))

        for gate in gate_objects['gates']:
            gate_type = gate['type']

            if gate_type == 'Polygon':
                gate['result'] = parse_polygon_gate(events, channel_labels, gate)
                gated_events = np.vstack((gated_events, gate['result']['gated_events']))
            elif gate_type == 'Boolean':
                # Need to save these as the referenced gates may not have been
                # calculated yet, we'll parse them after all the other ones
                boolean_gates.append(gate)
            else:
                raise ValueError("Unsupported gate type: %s" % gate_type)

        parse_gates(gated_events, channel_labels, gate_objects['children'])

    if len(boolean_gates) > 0:
        parse_boolean_gates(events, gating_dict, boolean_gates)


def parse_boolean_gates(events, gating_dict, boolean_gate_list):

    for gate in boolean_gate_list:
        specs = [g.strip() for g in gate['specification'].split('&')]
        group_inclusion = []
        for g in specs:
            if g[0] == '!':
                group_inclusion.append(False)
            else:
                group_inclusion.append(True)

        groups = []

        for g in gate['groups']:
            g_split = g.split('/')
            if len(g_split) > 2:
                warnings.warn("Nested boolean gates are not supported (%s)" % g)
                groups.append(None)
                continue

            groups.append(g_split[1])

        # start with including all events
        all_include_events = np.broadcast_to(True, events[:, 0].shape)

        for i, g in enumerate(groups):
            if group_inclusion[i] and g is not None:
                invert = False
            elif g is not None:
                invert = True  # invert because it's an exclusion gate
            else:
                continue

            for group_gate in gating_dict[g]['gates']:
                if len(group_gate['result']['gated_events']) > 0:
                    # get the event indices for this gate
                    gate_events = group_gate['result']['gated_events'][:, 0]

                    # find boolean array of these indices from all parent events
                    gate_include_events = np.in1d(
                        events[:, 0],
                        list(gate_events),
                        invert=invert
                    )
                    all_include_events = np.logical_and(
                        all_include_events,
                        gate_include_events
                    )

        gated_events = events[all_include_events]

        gate['result'] = {
            'gated_events': gated_events,
            'ungated_count': events.shape[0]
        }


def apply_gating_hierarchy(events, channel_labels, gating_dict):
    """
    Extract events from root gates and recurse on any children

    :param events: NumPy array of events on which to apply the gate
    :param channel_labels: dictionary of channel labels (keys are channel #'s)
    :param gating_dict: dictionary of gating hierarchy where current_gate
                        is a key at the root level. Gating dict is modified to
                        add the results to each gate.
    :return: None
    """

    parse_gates(events, channel_labels, gating_dict)
