import pandas as pd
from metaflow import FlowSpec, Parameter, step

from profiling.calculate_physchem import calculate_properties
from profiling.plot_physchem import plot_property_profiles


class Property_Profile_Generation(FlowSpec):
    library_name = Parameter(
        "library_name",
        help="The XLIB or XPL number of the library. For example, XLIB0001 or XPL0002",
        default="",
    )

    mode = Parameter(
        "mode", help='"presynthesis" or "postsynthesis"', default="postsynthesis"
    )

    @step
    def start(self):
        print(self.library_name)
        self.output_dir = f"/mnt/p/Discovery Chemistry/Library Synthesis/Enumeration/Library Enumeration/{self.library_name}/{self.mode}/output/"
        if self.mode == "postsynthesis":
            filename = self.output_dir + f"{self.library_name}_output.xlsx"
        else:
            filename = self.output_dir + f"{self.library_name}_DELflex_output.xlsx"
        self.df = pd.read_excel(
            filename,
            sheet_name="Instance Enumeration",
        )
        self.next(self.calculate_missing_properties)

    @step
    def calculate_missing_properties(self):
        self.df = calculate_properties(self.df)
        self.next(self.plot)

    @step
    def plot(self):
        plot_property_profiles(self.df, self.output_dir + "/property_profiles/")
        self.next(self.end)

    @step
    def end(self):
        pass


if __name__ == "__main__":
    Property_Profile_Generation()
