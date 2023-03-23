classdef superjog

   properties
      idx_seg_list
      type
      rn_list
      
   end

   methods
      function obj = MyClass(val)
         if nargin > 0
            obj.Prop = val;
         end
      end
   end
end