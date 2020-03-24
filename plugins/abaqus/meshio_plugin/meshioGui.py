"""GUI for meshio and related tools in Postprocessing."""
#  import time
from abaqusGui import (
    AFXSELECTFILE_ANY,
    DIALOG_ACTIONS_SEPARATOR,
    FRAME_GROOVE,
    FRAME_THICK,
    FXMAPFUNC,
    HSCROLLING_OFF,
    LIST_BROWSESELECT,
    LIST_MULTIPLESELECT,
    SEL_COMMAND,
    AFXComboBox,
    AFXDataDialog,
    AFXDialog,
    AFXFileSelectorDialog,
    AFXForm,
    AFXList,
    AFXStringTarget,
    AFXTextField,
    FXButton,
    FXCheckButton,
    FXGroupBox,
    FXHorizontalFrame,
    FXLabel,
    FXVerticalFrame,
    sendCommand,
    session,
    showAFXWarningDialog,
)

import re

class CreateSymmTensor(AFXDataDialog):
    """Create a new tensor field from given variables."""

    [ID_WARNING] = range(AFXForm.ID_LAST, AFXForm.ID_LAST + 1)

    def __init__(self, form):
        """Set up the initial dialog and connect Buttons to actions."""
        self.form = form
        # find current viewport
        currentViewport = session.viewports[session.currentViewportName]
        # assign odb file from current viewport
        self.odbFile = currentViewport.displayedObject
        # get file name and path
        self.odbFileNameFull = self.odbFile.path

        AFXDataDialog.__init__(
            self, form, "Create Tensor Field", 0, DIALOG_ACTIONS_SEPARATOR
        )

        FXLabel(self, "This plugin creates a tensor from state variables.")

        nvis = 0

        self.nameTextField = AFXTextField(self, 35, "Field name:")
        self.descTextField = AFXTextField(self, 35, "Description:")

        variables = []
        for var in currentViewport.odbDisplay.fieldVariables.variableList:
            variables.append(var[0])

        nvis = len(variables)
        self.s1combo = AFXComboBox(self, 0, nvis, "Component 11:")
        self.s2combo = AFXComboBox(self, 0, nvis, "Component 22:")
        self.s3combo = AFXComboBox(self, 0, nvis, "Component 33:")
        self.s4combo = AFXComboBox(self, 0, nvis, "Component 12:")
        self.s5combo = AFXComboBox(self, 0, nvis, "Component 23:")
        self.s6combo = AFXComboBox(self, 0, nvis, "Component 31:")
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

        self.appendActionButton("Create Field", self, self.ID_CLICKED_APPLY)
        self.appendActionButton(self.DISMISS)

        FXMAPFUNC(
            self, SEL_COMMAND, self.ID_CLICKED_APPLY, CreateSymmTensor.doCustomChecks
        )
        FXMAPFUNC(self, SEL_COMMAND, self.ID_WARNING, CreateSymmTensor.onCmdWarning)

    def doCustomChecks(self, sender, sel, ptr):
        """Ask user, because this is a permanent modification of the ODB."""
        showAFXWarningDialog(
            self,
            "This will re-load the ODB in write mode and"
            " adds fields permantenly. Are you sure?",
            AFXDialog.YES | AFXDialog.NO,
            self,
            self.ID_WARNING,
        )

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
        cmdstr = (
            "tensor('"
            + self.odbFileNameFull
            + "', '"
            + field_name
            + "', '"
            + description
            + "', '"
            + s1
            + "', '"
            + s2
            + "', '"
            + s3
            + "', '"
            + s4
            + "', '"
            + s5
            + "', '"
            + s6
            + "')"
        )
        sendCommand(cmdstr)
        self.form.deactivate()
        return 1


class CreateVector(AFXDataDialog):
    """Create a new vector field from given variables."""

    [ID_WARNING] = range(AFXForm.ID_LAST, AFXForm.ID_LAST + 1)

    def __init__(self, form):
        """Set up the initial dialog and connect Buttons to actions."""
        self.form = form
        # find current viewport
        currentViewport = session.viewports[session.currentViewportName]
        # assign odb file from current viewport
        self.odbFile = currentViewport.displayedObject
        # get file name and path
        self.odbFileNameFull = self.odbFile.path

        AFXDataDialog.__init__(
            self, form, "Create Vector Field", 0, DIALOG_ACTIONS_SEPARATOR
        )

        FXLabel(self, "This plugin creates a vector from state variables.")

        nvis = 0

        self.nameTextField = AFXTextField(self, 35, "Field name:")
        self.descTextField = AFXTextField(self, 35, "Description:")

        variables = []
        for var in currentViewport.odbDisplay.fieldVariables.variableList:
            variables.append(var[0])

        nvis = len(variables)
        self.s1combo = AFXComboBox(self, 0, nvis, "Component 1:")
        self.s2combo = AFXComboBox(self, 0, nvis, "Component 2:")
        self.s3combo = AFXComboBox(self, 0, nvis, "Component 3:")
        for variable in sorted(variables):
            self.s1combo.appendItem(variable)
            self.s2combo.appendItem(variable)
            self.s3combo.appendItem(variable)
        self.s1combo.setMaxVisible(10)
        self.s2combo.setMaxVisible(10)
        self.s3combo.setMaxVisible(10)

        self.appendActionButton("Create Field", self, self.ID_CLICKED_APPLY)
        self.appendActionButton(self.DISMISS)

        FXMAPFUNC(self, SEL_COMMAND, self.ID_CLICKED_APPLY, CreateVector.doCustomChecks)
        FXMAPFUNC(self, SEL_COMMAND, self.ID_WARNING, CreateVector.onCmdWarning)

    def doCustomChecks(self, sender, sel, ptr):
        """Ask user, because this is a permanent modification of the ODB."""
        showAFXWarningDialog(
            self,
            "This will re-load the ODB in write mode and"
            " adds fields permantenly. Are you sure?",
            AFXDialog.YES | AFXDialog.NO,
            self,
            self.ID_WARNING,
        )

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
        cmdstr = (
            "vector('"
            + self.odbFileNameFull
            + "', '"
            + field_name
            + "', '"
            + description
            + "', '"
            + s1
            + "', '"
            + s2
            + "', '"
            + s3
            + "')"
        )
        sendCommand(cmdstr)
        self.form.deactivate()
        return 1


class ExportODB(AFXDataDialog):
    """Export ODB with meshio."""

    [ID_CLICKED_FILE_BUTTON] = range(AFXDataDialog.ID_LAST, AFXDataDialog.ID_LAST + 1)

    def __init__(self, form):

        def atoi(text):
            return int(text) if text.isdigit() else text

        def natural_keys(text):
            '''
            alist.sort(key=natural_keys) sorts in human order
            http://nedbatchelder.com/blog/200712/human_sorting.html
            (See Toothy's implementation in the comments)
            '''
            return [ atoi(c) for c in re.split(r'(\d+)', text) ]
        
        """Set up the initial dialog and connect Buttons to actions."""
        self.form = form
        # find current viewport
        currentViewport = session.viewports[session.currentViewportName]
        # assign odb file from current viewport
        self.odb = currentViewport.displayedObject
        # get file name and path
        self.odbFileNameFull = self.odb.path

        self.file_name = AFXStringTarget()

        AFXDataDialog.__init__(self, form, "Export ODB", 0, DIALOG_ACTIONS_SEPARATOR)

        FXLabel(self, "This plugin exports the ODB with Meshio.")

        self.variables = []
        for var in currentViewport.odbDisplay.fieldVariables.variableList:
            self.variables.append(var[0])
        
        self.variables.sort(key=natural_keys)
        
        hf_selectors = FXHorizontalFrame(self)

        gb_instances = FXGroupBox(hf_selectors, "Instance", FRAME_GROOVE)
        gb_instances_label = FXVerticalFrame(gb_instances, FRAME_THICK)
        self.instlist = AFXList(
            gb_instances_label, 10, None, 0, LIST_BROWSESELECT | HSCROLLING_OFF
        )
        for instance in self.odb.rootAssembly.instances.keys():
            self.instlist.appendItem(instance)

        gb_frames = FXGroupBox(hf_selectors, "Frame", FRAME_GROOVE)
        gb_frames_labels = FXVerticalFrame(gb_frames, FRAME_THICK)
        self.framelist = AFXList(
            gb_frames_labels, 10, None, 0, LIST_BROWSESELECT | HSCROLLING_OFF
        )
        for step in currentViewport.odbDisplay.fieldSteps:
            step_name = step[0]
            frames = step[7]
            for i, f in enumerate(frames):
                self.framelist.appendItem("%s: %d (%s)" % (step_name, i, f))

        gb_variables = FXGroupBox(hf_selectors, "Variables", FRAME_GROOVE)
        gb_variables_label = FXVerticalFrame(gb_variables, FRAME_THICK)
        self.varlist = AFXList(
            gb_variables_label,
            10,
            None,
            0,
            LIST_BROWSESELECT | HSCROLLING_OFF | LIST_MULTIPLESELECT,
        )
        for i, var in enumerate(self.variables):
            self.varlist.appendItem(var, None, i)

        hf_file = FXHorizontalFrame(self)
        self.sourceTextField = AFXTextField(hf_file, 20, "Export to:")
        self.sourceButton = FXButton(
            hf_file, "Select File", None, self, self.ID_CLICKED_FILE_BUTTON
        )

        self.check_deform = FXCheckButton(self, "Export deformed geometry.")
        self.check_deform.setCheck(state=True)

        self.appendActionButton("Export", self, self.ID_CLICKED_APPLY)
        self.appendActionButton(self.DISMISS)

        FXMAPFUNC(self, SEL_COMMAND, self.ID_CLICKED_APPLY, ExportODB.export)
        FXMAPFUNC(self, SEL_COMMAND, self.ID_CLICKED_FILE_BUTTON, ExportODB.select_file)

    def select_file(self, sender, sel, ptr):
        """Create file dialog, when the coresponding button was hit."""
        # A pattern describes wich file types should be selected.
        patterns = (
            "VTK (*.vtk)\n"
            "VTK (*.vtu)\n"
            "STL (*.stl)\n"
            "Dolfin-XML (*.xml)\n"
            "Med (*.med)\n"
            "Medit (*.mesh)\n"
            "Gmsh4-binary (*.msh)\n"
            "Permas (*.post)\n"
            "Permas (*.post.gz)\n"
            "Permas (*.dato)\n"
            "Permas (*.dato.gz)\n"
            "Moab (*.h5m)\n"
            "Off (*.off)\n"
            "Xdmf (*.xdmf)\n"
            "Xmf (*.xmf)\n"
            "Mdpa (*.mdpa)\n"
            "SVG (*.svg)\n"
            "Patran (*.pat)\n"
            "Exodus (*.e)\n"
            "Exodus (*.ex2)\n"
            "Exodus (*.exo)\n"
            "All Files (*.*)"
        )

        # open a file dialog and set it's target
        self.fileDialog = AFXFileSelectorDialog(
            self,
            "Select Source File",
            self.file_name,
            None,
            AFXSELECTFILE_ANY,
            patterns,
        )
        self.fileDialog.create()
        self.fileDialog.show()
        return 1
    
    def processUpdates(self):
        """Update fields.

        This is triggered by UI after each update.
        """
        self.sourceTextField.setText(self.file_name.getValue())

    def export(self, sender, sel, ptr):
        """Process inputs."""
        instance_item = self.instlist.getSingleSelection()
        inst = self.instlist.getItemText(instance_item)

        frame_item = self.framelist.getSingleSelection()
        substring = self.framelist.getItemText(frame_item).split("(")[0]
        step, frame = substring.split(": ")

        deformed = bool(self.check_deform.getCheck())

        variables = []
        for i, var in enumerate(self.variables):
            if self.varlist.isItemSelected(i):
                variables.append(var)

        tgt = r"".join(self.file_name.getValue())

        sendCommand("import meshio")
        sendCommand("from abq_meshio.abq_meshio_converter " "import convertODBtoMeshio")
        sendCommand("odb = session.openOdb('%s')" % self.odbFileNameFull)
        sendCommand("instance = odb.rootAssembly.instances['%s']" % inst)
        sendCommand("frame = odb.steps['%s'].frames[%d]" % (step, int(frame)))
        sendCommand(
            "odb_mesh = convertODBtoMeshio(instance, frame,"
            " list_of_outputs=['%s'], deformed=%s)"
            % ("','".join(variables), str(deformed))
        )
        sendCommand("meshio.write('%s', odb_mesh, write_binary=False)" % tgt)
        self.form.deactivate()
        return 1
