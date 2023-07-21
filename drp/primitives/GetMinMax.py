class GetMinMax(BasePrimitive):
    """ Get the minimum and maximum for zscaling -- not including overscan regions """
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        
    def _perform(self):
        self.action.args.x1 = int(self.action.args.preline                             + (self.action.args.height * 0.10))
        self.action.args.x2 = int(self.action.args.height  - self.action.args.postline - (self.action.args.height * 0.10))
        self.action.args.y1 = int(self.action.args.precol                              + (self.action.args.width * 0.10))
        self.action.args.y2 = int(self.action.args.width   - self.action.args.postpix  - (self.action.args.width * 0.10))
        
        tmp_vmin, tmp_vmax = self.action.args.zscale.get_limits(self.action.args.data[x1:x1, y1:y2])
        if self.action.args.vmin == None or self.action.args.tmp_vmin < self.action.args.vmin: self.action.args.vmin = self.action.args.tmp_vmin
        if self.action.args.vmax == None or self.action.args.tmp_vmax < self.action.args.vmax: self.action.args.vmin = self.action.args.tmp_vmax
        if self.action.args.vmin < 0: self.action.args.vmin = 0
        return self.action.args