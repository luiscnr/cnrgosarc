from arc_mapinfo import ArcMapInfo


class ArcProcessing:
    def __init__(self, arc_opt, verbose):
        self.verbose = verbose
        self.ami = ArcMapInfo(arc_opt, False)
        self.width = self.ami.area_def.width
        self.height = self.ami.area_def.height

        section = 'PROCESSING'
        