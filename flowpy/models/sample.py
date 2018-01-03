import numpy as np
import os
import flowio
import flowutils


class Sample(object):
    """
    Represents an Flow Cytometry Standard (FCS) sample
    """
    def __init__(self, fcs_file_path, use_file_compensation=False, track_indices=False):
        """
        fcs_file_path: path to FCS file
        use_file_compensation: whether to use the compensation matrix within the FCS file
        
        Note: Retrieving events always gives the sub-sampled data. Use subsample_count=0 
        to analyze all the events.
        """

        self.flow_data = flowio.FlowData(fcs_file_path)
        self.event_count = self.flow_data.event_count
        try:
            self.acquisition_date = self.flow_data.text['date']
        except KeyError:
            self.acquisition_date = None

        self.fcs_path = fcs_file_path
        self.filename = os.path.basename(fcs_file_path)

        self.channels = {}
        for chan_num, labels in self.flow_data.channels.items():
            if track_indices:
                new_chan_num = str(int(chan_num) + 1)
            else:
                new_chan_num = chan_num

            self.channels[new_chan_num] = labels

        if track_indices:
            self.channels['1'] = {'PnN': 'Index'}

        # compensation & transform info
        # TODO: implement extraction of compensation from FCS metadata
        self.compensation = None
        self.fluoro_indices = None

        # sub-sample info
        self.subsample_indices = None
        self.subsample_count = None

        # retrieving event arrays, always based on sub-sampled array
        self.events_compensated = None
        self.events_transformed = None

        # convert events to NumPy array
        self.raw_events = np.reshape(
            self.flow_data.events,
            (-1, self.flow_data.channel_count)
        )
        if track_indices:
            self.raw_events = np.hstack(
                (np.array(range(self.flow_data.event_count))[:, np.newaxis], self.raw_events)
            )

    def get_channel_numbers_by_channel_labels(self):
        channel_map = {}

        for c, labels in self.channels.items():
            channel_map[labels['PnN']] = int(c)

        return channel_map

    def generate_subsample(self, subsample_count, random_seed, filter_neg_scatter=False):
        """
        Sub-samples FCS sample. Set subsample_count to 0 (zero) to use all events
        
        Retrieve sub-sampled indices via subsample_indices.
        Retrieve sub-sample events via get_subsample()
        """
        # use all events if subsample count is zero (or less)
        if subsample_count <= 0:
            subsample_count = self.event_count

        # If flagged, filter out events with negative scatter values.
        # To do that we need the channel annotations
        if filter_neg_scatter:
            scatter_indices = []
            for c, labels in self.channels.items():
                if labels['PnN'][0:3] in ['FSC', 'SSC']:
                    scatter_indices.append(int(c) - 1)

            is_neg = self.raw_events[:, scatter_indices] < 0
            is_neg = np.where(is_neg.any(True))[0]
        else:
            is_neg = []

        if self.event_count - len(is_neg) < subsample_count:
            # The sample has fewer events than requested,
            # but we'll go ahead and use what we've got
            # TODO: at least throw a warning when this happens or make a flag-able option?
            subsample_count = self.event_count - len(is_neg)

        # generate random indices for subsample
        # using a new RandomState with given seed
        shuffled_indices = np.arange(self.event_count)
        shuffled_indices = np.delete(shuffled_indices, is_neg)
        rng = np.random.RandomState()
        rng.seed(random_seed)
        rng.shuffle(shuffled_indices)

        # save indices
        self.subsample_indices = np.sort(shuffled_indices[:subsample_count])

    @property
    def events_subsampled(self):
        """
        Retrieve sub-sample FCS events using subsample_indices.
        
        subsample_indices can be set via a given random seed the generate_subsample() method.
        If generate_subsample() is not manually called prior to accessing subsample_events,
        then all events will be used (including negative scatter events).
        :return: NumPy array of sub-sampled events
        """
        if self.subsample_indices is None:
            self.generate_subsample(0, 1, filter_neg_scatter=False)

        subsample = self.raw_events[self.subsample_indices]

        return subsample

    def compensate_events(self, compensation_matrix):
        """
        Applies compensation matrix to the sub-sampled events
        
        The compensation matrix must be a comma delimited text file with
        a header row containing the corresponding PnN channel labels.

        Retrieve compensated events via subsample_compensated_events
        """
        # flowutils compensate() takes the plain matrix and indices as
        # separate arguments
        # (also note channel #'s vs indices)

        # TODO: if this is called more than once, need to reset downstream transforms to None

        fluoro_indices = []

        if compensation_matrix is not None:
            # remove any white space at the beginning or end
            compensation_matrix = compensation_matrix.strip()

            # split lines by \r or \n
            compensation_matrix = compensation_matrix.splitlines()

            # split values by \t
            compensation_matrix = [line.split('\t') for line in compensation_matrix]

            # convert header to channel indices
            header = compensation_matrix[0]
            compensation_matrix = compensation_matrix[1:]  # just the matrix
            compensation_matrix = np.array(compensation_matrix, dtype=np.float)

            channel_map = self.get_channel_numbers_by_channel_labels()

            try:
                for label in header:
                    fluoro_indices.append(channel_map[label] - 1)
            except KeyError:
                raise KeyError("Compensation matrix labels do not match FCS file!")
        else:
            # assume identity matrix
            non_fluoro_channels = [
                'FSC-A',
                'FSC-H',
                'FSC-W',
                'SSC-A',
                'SSC-H',
                'SSC-W',
                'Time',
                'Index'
            ]

            for c, labels in self.channels.items():
                if labels['PnN'] not in non_fluoro_channels:
                    fluoro_indices.append(int(c) - 1)

            fluoro_indices.sort()
            compensation_matrix = np.identity(len(fluoro_indices))

        comp_data = flowutils.compensate.compensate(
            self.events_subsampled,
            compensation_matrix,
            fluoro_indices
        )

        self.fluoro_indices = fluoro_indices
        self.events_compensated = comp_data

    def transform_logicle(self, t=262144, w=0.5):
        """
        The default pre-scale factor is 0.003

        Retrieve transformed data via events_transformed
        """
        # TODO: check t default, need to calculate dynamically?

        x_data = flowutils.transforms.logicle(
            self.events_compensated,
            self.fluoro_indices,
            t=t,
            w=w
        )

        self.events_transformed = x_data

    def transform_asinh(self, pre_scale=0.003):
        """
        Applies inverse hyperbolic sine transform on compensated data.

        The default pre-scale factor is 0.003

        Retrieve transformed data via events_transformed
        """
        x_data = flowutils.transforms.asinh(
            self.events_compensated,
            self.fluoro_indices,
            pre_scale=pre_scale
        )

        self.events_transformed = x_data
