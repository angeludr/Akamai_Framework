class GetBias(BasePrimitive):
    """ Get the bias from the overscan region """
    
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        
    def _perform(self):
        self.action.args.x1 = 0
        self.action.args.x2 = self.action.args.height
        self.action.args.y1 = self.action.args.width - self.action.args.postpix + 1
        self.action.args.y2 = self.action.args.width
        
        self.action.args.bias = np.median(self.action.args.data[x1:x2, y1:y2], axis=1)
        self.action.args.bias = np.array(self.action.args.bias, dtype=np.int64)
        return self.action.args