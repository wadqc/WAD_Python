import dicom


### Helper functions
def getValue(ds, label):
    """Return the value of a pydicom DataElement in Dataset identified by label.

    ds: pydicom Dataset
    label: dicom identifier, in either pydicom Tag object, string or tuple form.
    """
    if isinstance(label, str):
        try:
            #Assume form "0x0008,0x1030"
            tag = dicom.tag.Tag(label.split(','))
        except ValueError:
            try:
                #Assume form "SeriesDescription"
                tag = ds.data_element(label).tag
            except (AttributeError, KeyError):
                #`label` string doesn't represent an element of the DataSet
                return None
    else:
        #Assume label is of form (0x0008,0x1030) or is a pydicom Tag object.
        tag = dicom.tag.Tag(label)

    try:
        return str(ds[tag].value)
    except KeyError:
        #Tag doesn't exist in the DataSet
        return None


def isFiltered(ds, filters):
    """Return True if the Dataset `ds` complies to the `filters`,
    otherwise return False.
    """
    for tag, value in filters.items():
        if not str(getValue(ds, tag)) == str(value):
            # Convert both values to string before comparison. Reason is that
            # pydicom can return 'str', 'int' or 'dicom.valuerep' types of data.
            # Similarly, the user (or XML) supplied value can be of any type.
            return False
    return True


def flatten(nestedlist):
    """Flatten a list of lists"""
    return [item for sublist in nestedlist for item in sublist]


def applyFilters(series_filelist, filters):
    """Apply `filters` to the `series_filelist` and return the filtered list.

    First, convert `filters` from an ElementTree Element to a dictionary
    Next, create a new list in the same shape as `series_filelist`, but only
    include filenames for which isFiltered returns True.
    Only include sublists (i.e., series) which are non empty.
    """
    #Turn ElementTree element attributes and text into filters
    filter_dict = {element.attrib["name"]: element.text for element in filters}

    filtered_series_filelist = []
    #For each series in the series_filelist (or, study):
    for instance_filelist in series_filelist:
        #Filter filenames within each series
        filtered_instance_filelist = [fn for fn in instance_filelist
            if isFiltered(
                dicom.read_file(fn, stop_before_pixels=True), filter_dict)]
        #Only add the series which are not empty
        if filtered_instance_filelist:
            filtered_series_filelist.append(filtered_instance_filelist)

    return filtered_series_filelist


class PluginData(object):
    """Contain references to the dicom files available to the plugin and provide
    various methods to convieniently access the corresponding pydicom Dataset
    objects.

    Every 'get...' method can be called with keyword arguments which are passed
    to the dicom.read_file() function calls. This is useful, for instance, when
    only the headers of the files are needed by the plugin. E.g.:
    >> data.getSeriesByDescription("Flat Field", stop_before_pixels=True)
    """

    def __init__(self, series_filelist, filters):


        if filters:
            self.series_filelist = applyFilters(series_filelist, filters)
        else:
            #If no filters are provided, keep all files. This would happen with
            #the applyFilters function aswell, but that would unnecessarily load
            #all dicom files.
            self.series_filelist = series_filelist

    def loadDicomFunction(self, **kwargs):
        return lambda fn: dicom.read_file(fn, **kwargs)

    def getAllInstances(self, **kwargs):
        """Return a flattened list of all pydicom Dataset objects obtained from
        series_filelist
        """

        return flatten(self.getAllSeries(**kwargs))

    def getAllSeries(self, **kwargs):
        """Return a nested list of all pydicom Dataset objects obtained from
        series_filelist
        """
        func = self.loadDicomFunction(**kwargs)
        return [[func(f) for f in filelist] for filelist in self.series_filelist]

    def getSeriesByDescription(self, description, **kwargs):
        """Return the series for which the SeriesDescription matches
        `description` argument
        """
        func = self.loadDicomFunction(**kwargs)
        return [[func(f) for f in filelist] for filelist in self.series_filelist
            if description == dicom.read_file(filelist[0]).SeriesDescription]

    def getInstanceByTags(self, filters, **kwargs):
        """Return a list of dicom instances which satisfy the specified filters.

        filters: dictionary where each key, value pair is a pydicom DataElement
            name or tag number and the corresponding filter

        Example:
        myfilter = {
            "ExposureTime": "100",
            "PatientName": "WAD",
            (0x0018,0x1405): "13"
        }
        """
        func = self.loadDicomFunction(**kwargs)
        instances = (func(fn) for fn in flatten(self.series_filelist))
        return [ds for ds in instances if isFiltered(ds, filters)]



