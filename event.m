classdef event

   properties
      segtype
      rate
      idx_seg
      action
   end

   methods
      function obj = MyClass(val)
         if nargin > 0
            obj.Prop = val;
         end
      end
   end
end