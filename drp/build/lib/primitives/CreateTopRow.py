class CreateTopRow(BasePrimitive):
    """ Creating second row of detector mosaic -- ext#: 4 - 8 """
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        
    def _perform(self):
        self.action.args.top = []
        for self.action.args.arr in self.action.args.alldata[4:]:
            self.action.args.arr = np.flip(self.action.args.arr)
            self.action.args.top.append(self.action.args.arr)
        self.action.args.top = np.concatenate(self.action.args.top, axis=1)
        return self.action.args