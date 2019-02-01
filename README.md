### **Note: FlowPy is no longer being developed or maintained due to the difficulty in supporting multiple versions of the proprietary FlowJo workspace format and the lack of an XML schema definition for validating those documents. [FlowIO](https://github.com/whitews/flowio) & [FlowUtils](https://github.com/whitews/flowutils) are still being maintained for lower level interaction with FCS data, and a new library, [FlowKit](https://github.com/whitews/flowkit), is in development to support high-level interaction including full GatingML 2.0 support.**

# FlowPy

A Python library for high level interaction with flow cytometry data (FCS files), including support for FlowJo XML workspace files.

## Requirements

* matplotlib
* numpy
* pandas
* [flowio](https://github.com/whitews/FlowIO)
* [flowutils](https://github.com/whitews/FlowUtils)

## Usage

### Classes

FlowPy uses 2 basic classes for interacting with FCS data: a **Sample** class representing a single FCS file, and the **Workspace** class representing a FlowJo workspace (saved in XML format).

#### Sample class

A Sample class instance represents a single FCS file, providing several methods for common operations as well as attributes for commonly used metadata such as acquisition date and event count.

**Initialization**

`Sample(fcs_file_path, track_indices=False)`

Initialization of a Sample instance given the path to an FCS file. If True, `track_indices` adds an index column to the event data for tracking individual events over all analysis operations.

**Attributes**

`metadata`

Contains all the key / value pairs found in the given FCS file's metadata

`event_count`

The number of events found in the given FCS file

`fcs_path`

File path of given FCS file

`filename`

File name of the given FCS file (not the filename value found in the metadata)

`channels`

A dictionary containing the channel labels of the FCS file. The keys are channel numbers (not indices). The value is a dictionary where the required keyword 'PnN' is always present, the optional 'PnS' key is recorded, if present in the FCS file.

`subsample_indices`

A NumPy array containing the indices of the subsampled events. Note, this is only populated after the `generate_subsample` method is called, before which the value is `None`.

`subsample_count`

The number of events used for subsampling. Note, this is only populated after the `generate_subsample` method is called, before which the value is `None`.

`events_subsampled`

A property returning the sub-sampled events. If `generate_subsample` has not been called, the property will call the method to subsample with all events, a random seed of 1, and `filter_neg_scatter` set to False.

`events_compensated`

A NumPy array containing the compensated events. Note, this is only populated after the `compensate_events` method is called, before which the value is `None`.

`events_transformed`

A NumPy array containing the transformed events. Only compensated events are transformed, so if a transform is desired on uncompensated data, an identity matrix should be provided to the `compensate_events` method.

`raw_events`

A NumPy array containing the raw, unprocessed FCS event data. If the Sample instance was initialized using `track_indices=True`, then the raw events array will contain an index column as the first column.

**Methods**

`get_channel_numbers_by_channel_labels()`

Returns a dictionary mapping channel label keys to channel number values (not channel indices)

`generate_subsample(self, subsample_count, random_seed, filter_neg_scatter=False)`

Sub-samples event data in the FCS file by the given `subsample_count` using the given `random_seed`. The sub-sampled events are sorted by index after shuffling. If the requested subsample count is less than the number of events in the FCS file, the sub-sample will simply be all events. If `filter_neg_scatter` is True, any negative scatter events will be removed prior to sub-sampling.

Note: If `generate_subsample` has not been called prior to accessing the `events_subsampled` property, the subsample will be generated with all events, a random seed of 1, and `filter_neg_scatter` set to False.

`compensate_events(self, compensation_matrix)`

Applies `compensation_matrix` to the sub-sampled events. The compensation matrix must be a comma delimited text file with a header row containing the corresponding PnN channel labels.

Retrieve compensated events via `events_compensated` attribute

`transform_logicle(self, t=262144, w=0.5)`

Apply a logicle transform to the **compensated** events. Retrieve transformed data via `events_transformed`.

`transform_asinh(self, pre_scale=0.003)`

Apply a inverse hyperbolic sine transform to the **compensated** events. Retrieve transformed data via `events_transformed`.

#### Workspace class

A Sample class instance represents a single FCS file, providing several methods for common operations as well as attributes for commonly used metadata such as acquisition date and event count.

**Initialization**

`Workspace(xml_workspace)`

Initialization of a Workspace instance given the path to an FlowJo XML workspace file. The Workspace class contains gating strategies, of which there are two types (mirroring the structure of FlowJo workspaces): individual FCS sample strategies and those for groups of samples.

**Attributes**

`samples`

A dictionary containing information about individual FCS samples used in the workspace. The keys are the sample ID strings found in the XML workspace. The value is a dictionary with the following keys:

* `filename`: FCS file name found in XML workspace. Note: this value may differ from the actual filename on the filesystem as FlowJo defaults to using the `$FIL` metadata key value, which may or may not correspond to the actual filename if the original file has been renamed.
* `eventCount`: The total number of events found in the FCS file
* `populations`: A nested dictionary structure of the gating hierarchy. The key is the gate name, the value is a dictionary with keys `gates` and `children`. The `gates` value contains the gate boundaries and the `children` contain sub-populations with the same structure as the parent population.

`groups`

* `group_name`: Group name found in XML workspace for a collection of FCS samples.
* `samples`: The sample IDs of the samples within the group
* `populations`: A nested dictionary structure of the gating hierarchy. The key is the gate name, the value is a dictionary with keys `gates` and `children`. The `gates` value contains the gate boundaries and the `children` contain sub-populations with the same structure as the parent population.

**Methods**

`get_gate_hierarchies()`

Returns a dictionary listing the sample and group gate hierarchy names by ID.

`find_matching_gate_hierarchies(fcs_file_path)`

Returns a dictionary of gate hierarchies found within the workspace to which the given FCS file belongs

`analyze_sample(self, fcs_file_path, comp_matrix, gate_type, gate_id)`

Applies a gating hierarchy to the given FCS file. Returns results as a dictionary with the following keys:

* `filename`: The file name of the given FCS file
* `populations`: A nested dictionary structure of the gating hierarchy similar to the populations dictionary of the `samples` and `groups` attributes, except with the addition of the event counts from the FCS file calculated for each gate.
* `report`: A Pandas DataFrame containing a summary report of the gate populations, with the columns 'parent_path', 'label', 'type', 'parent_count', 'count', and 'relative_percent'.
