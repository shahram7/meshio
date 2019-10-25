"""GUI for meshio and related tools in Postprocessing."""

from abaqusGui import (AFXDataDialog, AFXComboBox, session, FXLabel, FXMAPFUNC,
                       sendCommand, AFXTextField, DIALOG_ACTIONS_SEPARATOR,
                       SEL_COMMAND, showAFXWarningDialog, AFXDialog, AFXForm)
from abaqusConstants import SCALAR


class CreateSymmTensor(AFXDataDialog):
    """Create a new tensor field from given variables."""

    [
        ID_WARNING,
    ] = range(AFXForm.ID_LAST, AFXForm.ID_LAST+1)

    def __init__(self, form):
        """Set up the initial dialog and connect Buttons to actions."""
        self.form = form
        # find current viewport
        currentViewport = session.viewports[session.currentViewportName]
        # assign odb file from current viewport
        self.odbFile = currentViewport.displayedObject
        # get file name and path
        self.odbFileNameFull = self.odbFile.path

        AFXDataDialog.__init__(self,
                               form,
                               'Create Tensor Field',
                               0,
                               DIALOG_ACTIONS_SEPARATOR)

        FXLabel(self, 'This plugin creates a tensor from state variables.')

        nvis = 0

        self.nameTextField = AFXTextField(self, 35, 'Field name:')
        self.descTextField = AFXTextField(self, 35, 'Description:')

        variables = set()
        step_name = self.odbFile.steps.keys()[0]
        frame = self.odbFile.steps[step_name].frames[0]
        for name in frame.fieldOutputs.keys():
            if frame.fieldOutputs[name].type == SCALAR:
                variables.add(name)
        variables = list(variables)

        nvis = len(variables)
        self.s1combo = AFXComboBox(self, 0, nvis, 'Component 11:')
        self.s2combo = AFXComboBox(self, 0, nvis, 'Component 22:')
        self.s3combo = AFXComboBox(self, 0, nvis, 'Component 33:')
        self.s4combo = AFXComboBox(self, 0, nvis, 'Component 12:')
        self.s5combo = AFXComboBox(self, 0, nvis, 'Component 23:')
        self.s6combo = AFXComboBox(self, 0, nvis, 'Component 31:')
        for variable in sorted(variables):
            self.s1combo.appendItem(variable)
            self.s2combo.appendItem(variable)
            self.s3combo.appendItem(variable)
            self.s4combo.appendItem(variable)
            self.s5combo.appendItem(variable)
            self.s6combo.appendItem(variable)
        self.s1combo.setMaxVisible(10)
        self.s2combo.setMaxVisible(10)
        self.s3combo.setMaxVisible(10)
        self.s4combo.setMaxVisible(10)
        self.s5combo.setMaxVisible(10)
        self.s6combo.setMaxVisible(10)

        self.appendActionButton('Create Field', self, self.ID_CLICKED_APPLY)
        self.appendActionButton(self.DISMISS)

        FXMAPFUNC(self, SEL_COMMAND, self.ID_CLICKED_APPLY,
                  CreateSymmTensor.doCustomChecks)
        FXMAPFUNC(self, SEL_COMMAND, self.ID_WARNING,
                  CreateSymmTensor.onCmdWarning)

    def doCustomChecks(self, sender, sel, ptr):
        """Ask user, because this is a permanent modification of the ODB."""
        showAFXWarningDialog(self,
                             'This will re-load the ODB in write mode and'
                             ' adds fields permantenly. Are you sure?',
                             AFXDialog.YES | AFXDialog.NO,
                             self, self.ID_WARNING)

    def onCmdWarning(self, sender, sel, ptr):
        """Send command to build tensor or abort."""
        if sender.getPressedButtonId() == AFXDialog.ID_CLICKED_YES:
                self.create_field()
        elif sender.getPressedButtonId() == AFXDialog.ID_CLICKED_NO:
                self.form.deactivate()

    def create_field(self):
        """Process inputs."""
        field_name = self.nameTextField.getText()
        description = self.descTextField.getText()

        item = self.s1combo.getCurrentItem()
        s1 = self.s1combo.getItemText(item)
        item = self.s2combo.getCurrentItem()
        s2 = self.s2combo.getItemText(item)
        item = self.s3combo.getCurrentItem()
        s3 = self.s3combo.getItemText(item)
        item = self.s4combo.getCurrentItem()
        s4 = self.s4combo.getItemText(item)
        item = self.s5combo.getCurrentItem()
        s5 = self.s5combo.getItemText(item)
        item = self.s6combo.getCurrentItem()
        s6 = self.s6combo.getItemText(item)

        sendCommand("from abq_meshio.combine import *")
        cmdstr = ("tensor('" + self.odbFileNameFull + "', '" +
                  field_name + "', '" +
                  description + "', '" +
                  s1 + "', '" +
                  s2 + "', '" +
                  s3 + "', '" +
                  s4 + "', '" +
                  s5 + "', '" +
                  s6 + "')")
        sendCommand(cmdstr)
        self.form.deactivate()
        return 1


class CreateVector(AFXDataDialog):
    """Create a new tensor field from given variables."""

    [
        ID_WARNING,
    ] = range(AFXForm.ID_LAST, AFXForm.ID_LAST+1)

    def __init__(self, form):
        """Set up the initial dialog and connect Buttons to actions."""
        self.form = form
        # find current viewport
        currentViewport = session.viewports[session.currentViewportName]
        # assign odb file from current viewport
        self.odbFile = currentViewport.displayedObject
        # get file name and path
        self.odbFileNameFull = self.odbFile.path

        AFXDataDialog.__init__(self,
                               form,
                               'Create Vector Field',
                               0,
                               DIALOG_ACTIONS_SEPARATOR)

        FXLabel(self, 'This plugin creates a vector from state variables.')

        nvis = 0

        self.nameTextField = AFXTextField(self, 35, 'Field name:')
        self.descTextField = AFXTextField(self, 35, 'Description:')

        variables = set()
        step_name = self.odbFile.steps.keys()[0]
        frame = self.odbFile.steps[step_name].frames[0]
        for name in frame.fieldOutputs.keys():
            if frame.fieldOutputs[name].type == SCALAR:
                variables.add(name)
        variables = list(variables)

        nvis = len(variables)
        self.s1combo = AFXComboBox(self, 0, nvis, 'Component 1:')
        self.s2combo = AFXComboBox(self, 0, nvis, 'Component 2:')
        self.s3combo = AFXComboBox(self, 0, nvis, 'Component 3:')
        for variable in sorted(variables):
            self.s1combo.appendItem(variable)
            self.s2combo.appendItem(variable)
            self.s3combo.appendItem(variable)
        self.s1combo.setMaxVisible(10)
        self.s2combo.setMaxVisible(10)
        self.s3combo.setMaxVisible(10)

        self.appendActionButton('Create Field', self, self.ID_CLICKED_APPLY)
        self.appendActionButton(self.DISMISS)

        FXMAPFUNC(self, SEL_COMMAND, self.ID_CLICKED_APPLY,
                  CreateVector.doCustomChecks)
        FXMAPFUNC(self, SEL_COMMAND, self.ID_WARNING,
                  CreateVector.onCmdWarning)

    def doCustomChecks(self, sender, sel, ptr):
        """Ask user, because this is a permanent modification of the ODB."""
        showAFXWarningDialog(self,
                             'This will re-load the ODB in write mode and'
                             ' adds fields permantenly. Are you sure?',
                             AFXDialog.YES | AFXDialog.NO,
                             self, self.ID_WARNING)

    def onCmdWarning(self, sender, sel, ptr):
        """Send command to build tensor or abort."""
        if sender.getPressedButtonId() == AFXDialog.ID_CLICKED_YES:
                self.create_field()
        elif sender.getPressedButtonId() == AFXDialog.ID_CLICKED_NO:
                self.form.deactivate()

    def create_field(self):
        """Process inputs."""
        field_name = self.nameTextField.getText()
        description = self.descTextField.getText()

        item = self.s1combo.getCurrentItem()
        s1 = self.s1combo.getItemText(item)
        item = self.s2combo.getCurrentItem()
        s2 = self.s2combo.getItemText(item)
        item = self.s3combo.getCurrentItem()
        s3 = self.s3combo.getItemText(item)

        sendCommand("from abq_meshio.combine import *")
        cmdstr = ("vector('" + self.odbFileNameFull + "', '" +
                  field_name + "', '" +
                  description + "', '" +
                  s1 + "', '" +
                  s2 + "', '" +
                  s3 + "')")
        sendCommand(cmdstr)
        self.form.deactivate()
        return 1
