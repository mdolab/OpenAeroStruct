import openmdao.api as om


class MultiCD(om.ExplicitComponent):
    def initialize(self):
        self.options.declare("n_points", types=int)

    def setup(self):
        self.n_points = self.options["n_points"]
        for i in range(self.n_points):
            self.add_input(str(i) + "_CD", val=0.0)

        self.add_output("CD", val=0.0)
        self.declare_partials("*", "*", val=1.0)

    def compute(self, inputs, outputs):
        outputs["CD"] = 0.0
        for i in range(self.n_points):
            outputs["CD"] += inputs[str(i) + "_CD"]
