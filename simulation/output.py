class Output:
    def __init__(
        self,
        output_folder: str = None,                 # Folder Output
        paraxial_model: bool = None,
        drs: bool = None
    ) -> None:
        if output_folder == None:
            self.output_folder = None
            self.SAVE = False
        else:
            self.output_folder = output_folder
            self.SAVE = True

        if paraxial_model == None:
            self.paraxial_model = False
        else:
            self.paraxial_model = paraxial_model

        if drs == None:
            self.drs = False
        else:
            self.drs = drs

    def get_output_folder(self):
        return self.output_folder
