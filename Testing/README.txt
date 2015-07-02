To add a new plugin to the PluginTesting framework, the following is required:

  1. The python sources of the new plugin (say NewPlugin) are in its own subfolder of plugins: e.g. Plugins/ModalityX/NewPlugin
  2. There is a corresponding Testing folder (e.g. Testing/ModalityX/NewPlugin) containing:
    a. an up-to-date config file 'config_something.xml' (e.g. Testing/ModalityX/NewPlugin/config_uniformity.xml)
    b. optionally a corresponding inputfile 'input_something.xml' (e.g. Testing/ModalityX/NewPlugin/input_uniformity.xml), which contains at least <analyselevel>level</analyselevel> where 'level' is series or study or instance; if this file is missing <analyselevel>series</analyselevel> is assumed
  3. There is a corresponding dicom data folder in the 'dropbox share' (e.g. Testing/ModalityX/NewPlugin/dicom_uniformity/)
  
Note:
  * the folder structures for Plugins and Testing and dicom data should be consistent
  * the identifier of the test should be consistent in the naming of the config file, and the input file and the dicom data folder

