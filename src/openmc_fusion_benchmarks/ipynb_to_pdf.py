import os


class Report:
    def __init__(self):
        """Initializes the Report class with a list of files to convert."""
        self.notebook_files = []

    def add_ipynb_file(self, notebook_path: str):
        """Adds the path of a jupyter notebook file to the list."""
        if notebook_path.endswith(".ipynb"):
            self.notebook_files.append(notebook_path)
            print(f"File added: {notebook_path}")
        else:
            print(f"Error: {notebook_path} is not a .ipynb file.")

    def convert_to_pdf(self):
        """Converts all jupyter notebook files in the list to PDF."""
        if not self.notebook_files:
            print("Error: No files to convert.")
            return

        for notebook_path in self.notebook_files:
            print(f"Converting  {notebook_path} to PDF...")
            os.system(f"jupyter nbconvert --to pdf --no-input {notebook_path}")
        print("Conversion completed.")

    def modify(self):
        """Method for future modification. Raises a NotImplementedError."""
        raise NotImplementedError(
            "The modify method is not implemented yet and will be added in the future."
        )
