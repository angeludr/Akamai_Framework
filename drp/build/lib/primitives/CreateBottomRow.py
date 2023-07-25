class CreateBottomRow(BasePrimitive):
    """ Creating first row of detector mosaic -- ext#: 1 - 4 """
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        
    def _pre_condition(self):
        
        if self.action.args.handled < 8: 
            return False
        else:
            return True
    
    def _perform(self):
        self.action.args.bottom = np.concatenate(self.action.args.alldata[:4], axis=1)
        self.action.args.bottom = np.flipud(self.action.args.bottom)
        return self.action.args