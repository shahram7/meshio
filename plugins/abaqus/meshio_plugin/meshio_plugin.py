"""Provide tools to export results with meshio."""
import meshioGui as gui
from abaqusGui import AFXForm, getAFXApp

# Get toolset from Abaqus GUI application. The new plugins must be registered
# to this toolset in order to load at Abaqus startup
mainWindow = getAFXApp().getAFXMainWindow()
toolset = mainWindow.getPluginToolset()

# For each Plugin there is a new class derivated from AFXForm. Each constructor
# must call the base class constructor in order to pass the owner. Furthermore
# each class has a overridden function getFirstDialog, which defines the Dialog
# to start from this plugin. All these dialogs are implemented in
# fibermapGui.py


class CreateSymmTensor(AFXForm):
    """Form to build a new Tensor from existing state variables."""

    def __init__(self, owner):
        """Pass the toolset for initialization."""
        AFXForm.__init__(self, owner)

    def getFirstDialog(self):
        """Fire up the first dialog."""
        return gui.CreateSymmTensor(self)


class CreateVector(AFXForm):
    """Form to build a new Vector from existing state variables."""

    def __init__(self, owner):
        """Pass the toolset for initialization."""
        AFXForm.__init__(self, owner)

    def getFirstDialog(self):
        """Fire up the first dialog."""
        return gui.CreateVector(self)


class ExportODB(AFXForm):
    """Form to export ODB with meshio."""

    def __init__(self, owner):
        """Pass the toolset for initialization."""
        AFXForm.__init__(self, owner)

    def getFirstDialog(self):
        """Fire up the first dialog."""
        return gui.ExportODB(self)


toolset.registerGuiMenuButton(
    buttonText="Meshio|Create new symmetric tensor",
    object=CreateSymmTensor(toolset),
    version="1.0",
    author="Nils Meyer",
)

toolset.registerGuiMenuButton(
    buttonText="Meshio|Create new vector",
    object=CreateVector(toolset),
    version="1.0",
    author="Nils Meyer",
)

toolset.registerGuiMenuButton(
    buttonText="Meshio|Export ODB",
    object=ExportODB(toolset),
    version="1.0",
    author="Nils Meyer",
)
