import os
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET


def parseConfig(configfile):
    """Return a list of dicts, where each dict represents an action as specified
    in the config XML.

    Find every 'action' element in the ElementTree object. For each action, find
    the plugin, function and default_level (optional) text. If Elements for
    filters, limits and/or params are available, add the corresponding
    ElementTree Element object, otherwise add an empty list (to avoid None type)
    """
    xmlroot = ET.parse(configfile).getroot()

    action_list = []
    for action in xmlroot.findall('action'):
        entry = {}
        entry['plugin'] = action.find('plugin').text
        entry['function'] = action.find('function').text
        try:
            entry['default_level'] = action.find('default_level').text
        except AttributeError:
            entry['default_level'] = None

        entry['filters'] = action.find('filters') or []
        entry['limits'] = action.find('limits') or []
        entry['params'] = action.find('params') or []

        action_list.append(entry)

    return action_list


def parseInput(inputfile):
    """Return a dict representing several entries from the input XML (generated
    by the WAD server). Convert relative to absolute paths.

    Assume one patient and one study. Find all series and for each series add
    a list of its instance filepaths
    """
    xmlroot = ET.parse(inputfile).getroot()

    parsed_input = {}
    parsed_input['analysemodule_cfg'] = os.path.abspath(
        xmlroot.find('analysemodule_cfg').text)
    parsed_input['analysemodule_output'] = os.path.abspath(
        xmlroot.find('analysemodule_output').text)
    parsed_input['analyselevel'] = xmlroot.find('analyselevel').text
    parsed_input['serieslist'] = []

    for series in xmlroot.find('patient').find('study').findall('series'):
        instancelist = [instance.find('filename').text
            for instance in series.findall('instance')]
        parsed_input['serieslist'].append(
            [os.path.abspath(fn) for fn in instancelist]
        )

    return parsed_input

